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
 * Implementation of CAD-based mesh deformation.
 */
#pragma once

#include <vector>
#include <cmath>
#include <memory>
#include <t8_types/t8_vec.hxx>

/** 
 * The available RBF function types. 
 */
typedef enum { T8_RBF_CP_C2 = 0, T8_RBF_TPS, T8_RBF_FUNCTION_COUNT } t8_rbf_function_type;

struct t8_rbf_node
{
  /**The global ID of the node. */
  t8_gloidx_t global_id;
  /** The old position of the node. */
  t8_3D_vec position;
  /** The displacement. 
   * For boundary nodes this is known and calculated from the different coordinates of the CAD geometries. 
   * For inner nodes this is the result of the interpolation. */
  t8_3D_vec displacement;
};

/**
 * This struct will be a single boundary node from the mesh deformation. 
 * The new position of this node can be directly calculated from the new incoming CAD geometry given. 
 */
struct t8_rbf_boundary_node: public t8_rbf_node
{
  /* The calculated RBF coefficient alpha. This is needed for the interpolation to distribute the influence of the boundary nodes. */
  t8_3D_vec weight;
};

struct t8_rbf_function
{
  /** Destructor. */
  virtual ~t8_rbf_function () {};
  virtual double
  evaluate (double distance) const
    = 0;
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
  /** Constructor. 
   * \param[in] r   The support radius. Beyond the radius the function is zero and so there is no impact from this node movement.
   */
  t8_rbf_cpc2 (double r): radius (r)
  {
  }
  /** Destructor. */
  ~t8_rbf_cpc2 () {};
  /**
   * Solve the radial basis function.
   * Formula: psi(x) = (1 - x)^4 * (4x + 1) for (1-x) > 0.
   * where x = distance / support radius. The distance is the euclidean distance between two points. 
   * \param[in] distance  The euclidean distance between two points.
   * \return The function value psi. It returns 0.0 if the node is out of the chosen radius and so has no impact.   
   */
  double
  evaluate (double distance) const override
  {
    double r = distance / radius;
    if (r < 1.0) {
      double d = 1.0 - r;
      return (d * d * d * d) * (4.0 * r + 1.0);
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

 private:
  /** The chosen support radius. */
  double radius;
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
   * Solve the radial basis function.
   * Formula: psi(x) = x^2 * log(x).
   * where x is the euclidean distance between two points. 
   * \param[in] distance  The euclidean distance between two points.
   * \return The function value psi.  
   */
  double
  evaluate (double distance) const override
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
struct t8_rbf
{
 public:
  /** Constructor. 
   * \param[in] support_radius  The radius of the compact local support for the Wendland CP C2 function.
   *                            If the TPS RBF is used, we will not use the support_radius due to its global support.
   * \param[in] rbf_type        The RBF type to be used.
  */
  t8_rbf (double support_radius, t8_rbf_function_type rbf_type)
  {

    if (rbf_type == T8_RBF_CP_C2) {
      rbf_function = std::make_unique<t8_rbf_cpc2> (support_radius);
    }
    else if (rbf_type == T8_RBF_TPS) {
      rbf_function = std::make_unique<t8_rbf_tps> ();
    }
    else {
      t8_errorf ("ERROR: RBF attribute missing or not correct\n.");
      SC_ABORTF ("Unsupported RBF type.");
    }
  };

  /** Destructor. */
  ~t8_rbf () {};

  void
  set_boundary_nodes (std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node>&& boundary_node_data)
  {
    boundary_nodes.clear ();
    boundary_nodes.reserve (boundary_node_data.size ());
    for (auto& [global_vertex_id, boundary_node] : boundary_node_data) {
      boundary_nodes.push_back (std::move (boundary_node));
    }
  }

  /**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve ()
  {
  }

  /**
   * Interpolates the inner node.
   * \param[in, out] inner_node   The inner node which will be interpolated.
   */
  void
  interpolate (t8_rbf_node& inner_node) const
  {
    /** Reset the inner_node_displacement to zero. */
    inner_node.displacement.fill (0.0);
    /** Iterate over all boundary nodes. */
    for (const auto& boundary_node : boundary_nodes) {

      double distance = t8_dist (inner_node.position, boundary_node.position);
      const double psi = rbf_function->evaluate (distance);
      /** Check if the basis function value (the infleunce factor) is not equal to zero with a numerical tolerance of 1e-12. */
      if (std::abs (psi) > 1e-12) {

        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          inner_node.displacement[coordinate] += boundary_node.weight[coordinate] * psi;
        }
      }
    }
  }

 private:
  /** List of all registered support nodes. */
  std::vector<t8_rbf_boundary_node> boundary_nodes;
  /** */
  std::unique_ptr<t8_rbf_function> rbf_function;
};
