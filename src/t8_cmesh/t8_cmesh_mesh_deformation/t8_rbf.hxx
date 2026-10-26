/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_rbf.hxx
 * Radial basis function (RBF) used to interpolate the boundary displacements of a mesh deformation to the inner nodes.
 */
#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>
#include <t8_vector/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_cmesh_mesh_deformation_method.hxx>

/** Support radius used for globally supported RBFs, which have no finite support. */
constexpr double T8_RBF_GLOBAL_SUPPORT = std::numeric_limits<double>::infinity ();

/**
 * The available RBF function types.
 */
typedef enum {
  T8_RBF_CPC2 = 0,      /**< Compactly supported radial basis function (CPC2). */
  T8_RBF_TPS,           /**< Globally supported radial basis function (TPS). */
  T8_RBF_FUNCTION_COUNT /**< The number of available RBF function types. */
} t8_rbf_function_type;

/**
 * This struct will be a single boundary node from the mesh deformation.
 * The new position of this node can be directly calculated from the new incoming CAD geometry given.
 */
struct t8_rbf_boundary_node: public t8_mesh_deformation_node
{
  /** The calculated RBF coefficient alpha. This is needed for the interpolation to distribute the influence of the boundary nodes.
  * This weight represents the fixed influence assigned to each boundary node. */
  t8_3D_vec weight;
  /** The local support radius for the  boundary node. Only used for compactly supported RBFs. */
  double local_support_radius = T8_RBF_GLOBAL_SUPPORT;
};

/** A base struct representing a radial basis function. */
struct t8_rbf_function
{
  /** Destructor. */
  virtual ~t8_rbf_function () {};

  /**
   * Evaluate the radial basis function.
   * \param[in] distance The euclidean distance between two points.
   * \param[in] radius The support radius. Defaults to T8_RBF_GLOBAL_SUPPORT for globally supported RBFs.
   * \return The function value.
   */
  virtual double
  evaluate (const double distance, const double radius = T8_RBF_GLOBAL_SUPPORT) const
    = 0;

  /**
   * Check if the radial basis function is compactly supported.
   * \return true if the function is compactly supported, false otherwise.
   */
  virtual bool
  is_compactly_supported () const
    = 0;
};

/**
 * The CP C2 radial basis function which is compactly supported and will be used in the
 * mesh deformation to move the inner nodes from the known movement of the boundary nodes.
 */
struct t8_rbf_cpc2: public t8_rbf_function
{
  /** Constructor. */
  t8_rbf_cpc2 ()
  {
  }
  /** Destructor. */
  ~t8_rbf_cpc2 () {};
  /**
   * Evaluate the radial basis function.
   * Formula: psi(x) = (1 - x)^4 * (4x + 1) for (1-x) > 0,
   * where x = distance / support radius. The distance is the euclidean distance between two points.
   * \param [in] distance  The euclidean distance between two points.
   * \param [in] radius The support radius.
   * \return The function value psi. It returns 0.0 if the node is out of the chosen radius and so has no impact.
   */
  double
  evaluate (const double distance, const double radius) const override
  {
    /* CPC2 is compactly supported and requires a finite, positive support radius. */
    T8_ASSERT (std::isfinite (radius) && radius > 0.0);
    const double normalized_distance = distance / radius;
    if (normalized_distance < 1.0) {
      return std::pow (1.0 - normalized_distance, 4) * (4.0 * normalized_distance + 1.0);
    }
    return 0.0;
  }

  /**
   * The CPC2 is a compactly supported RBF.
   * \return true.
   */
  bool
  is_compactly_supported () const override
  {
    return true;
  }
};

/**
 * The TPS radial basis function with global support.
 */
struct t8_rbf_tps: public t8_rbf_function
{
  /** Constructor. The TPS RBF does not have a support radius, because it is globally supported.
 */
  t8_rbf_tps ()
  {
  }
  /** Destructor. */
  ~t8_rbf_tps () {};
  /**
   * Evaluate the radial basis function.
   * Formula: psi(x) = x^2 * log(x),
   * where x is the euclidean distance between two points.
   * \param[in] distance  The euclidean distance between two points.
   * \param [in] radius The support radius (not used for TPS because it is globally supported).
   * \return The function value psi.
   */
  double
  evaluate (const double distance, [[maybe_unused]] const double radius) const override
  {
    if (distance < 1e-12) {
      return 0.0;
    }
    return distance * distance * std::log (distance);
  }

  /**
   * The TPS is a globally supported RBF.
   * \return false.
   */
  bool
  is_compactly_supported () const override
  {
    return false;
  }
};

/**
 * Struct for mesh deformation using Radial Basis Functions (RBF).
 * It handles the interpolation of the boundary displacements to the inner nodes.
 */
struct t8_rbf: public t8_mesh_deformation_method
{
 public:
  /** Constructor.
   * \param [in] rbf_type        The mesh deformation/ RBF type to be used.
   * \param [in] scale_factor_support_radius      The scaling factor for the support radius of the CPC2 deformation. Note: Default value is 1.5 so it only has to be specified for CPC2 use.
  */
  t8_rbf (const t8_rbf_function_type rbf_type, const double scale_factor_support_radius = 1.5)
    : scale_factor_support_radius (scale_factor_support_radius)
  {

    if (rbf_type == T8_RBF_CPC2) {
      rbf_function = std::make_unique<t8_rbf_cpc2> ();
    }
    else if (rbf_type == T8_RBF_TPS) {
      rbf_function = std::make_unique<t8_rbf_tps> ();
    }
    else {
      SC_ABORT ("ERROR: RBF attribute missing or not correct. Unsupported RBF type.\n");
    }
  };

  /** Destructor. */
  ~t8_rbf () {};

  /**
   * Transfers the boundary node data to the internal data structure of the RBF.
   * \param [in] boundary_node_data Map of global vertex IDs to boundary nodes.
   */
  void
  set_boundary_nodes (const std::unordered_map<t8_gloidx_t, t8_mesh_deformation_node>& boundary_node_data) override
  {
    boundary_nodes.clear ();
    boundary_nodes.reserve (boundary_node_data.size ());
    for (const auto& [global_vertex_id, boundary_node] : boundary_node_data) {

      t8_rbf_boundary_node rbf_boundary_node;
      static_cast<t8_mesh_deformation_node&> (rbf_boundary_node) = boundary_node;
      rbf_boundary_node.global_id = global_vertex_id;
      rbf_boundary_node.weight.fill (0.0);

      /* Calculate the local support radius based on the mean neighbor distance. Globally supported RBFs keep the default T8_RBF_GLOBAL_SUPPORT. */
      if (rbf_function->is_compactly_supported ()) {
        rbf_boundary_node.local_support_radius = scale_factor_support_radius * rbf_boundary_node.mean_neighbor_distance;
      }

      boundary_nodes.push_back (std::move (rbf_boundary_node));
    }
  }

  /**
   * Solves the linear system A * alpha = displacements to find the weight.
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve () override;

  /**
   * Interpolates the inner node.
   * \param [in, out] inner_node   The inner node which will be interpolated.
   */
  void
  interpolate (t8_mesh_deformation_node& inner_node) const override
  {
    /* Reset the inner_node_displacement to zero. */
    inner_node.displacement.fill (0.0);
    /* Iterate over all boundary nodes. */
    for (const auto& boundary_node : boundary_nodes) {

      const double distance = t8_dist (inner_node.position, boundary_node.position);

      const double psi = rbf_function->evaluate (distance, boundary_node.local_support_radius);
      /* Check if the basis function value (the influence factor) is not equal to zero with a numerical tolerance of 1e-12. */
      if (std::abs (psi) > 1e-12) {

        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          /* Update the displacement for each coordinate with the weighted influence of the boundary node. */
          inner_node.displacement[coordinate] += boundary_node.weight[coordinate] * psi;
        }
      }
    }
  }

 private:
  /** List of all registered support nodes. */
  std::vector<t8_rbf_boundary_node> boundary_nodes;
  /** Pointer to the radial basis function. */
  std::unique_ptr<t8_rbf_function> rbf_function;
  /** The scaling factor for the support radius of compactly supported RBFs. */
  double scale_factor_support_radius;
};
