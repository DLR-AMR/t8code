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
#include <unordered_map>
#include <t8_types/t8_vec.hxx>
#include <Eigen/Sparse>
#include <Eigen/Dense>

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
  /** The calculated RBF coefficient alpha. This is needed for the interpolation to distribute the influence of the boundary nodes. 
  * This weight represents the fixed influence assigned to each boundary node. */
  t8_3D_vec weight;
  /** The local support radius for the  boundary node. */
  double local_support_radius;
};

struct t8_rbf_function
{
  /** Destructor. */
  virtual ~t8_rbf_function () {};
  virtual double
  evaluate (double distance, double radius) const
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
  /** Constructor. */
  t8_rbf_cpc2 ()
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
  evaluate (double distance, double radius) const override
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
  evaluate (double distance, [[maybe_unused]] double radius) const override
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
   * \param[in] rbf_type        The RBF type to be used.
  */
  t8_rbf (t8_rbf_function_type rbf_type)
  {

    if (rbf_type == T8_RBF_CP_C2) {
      rbf_function = std::make_unique<t8_rbf_cpc2> ();
    }
    else if (rbf_type == T8_RBF_TPS) {
      rbf_function = std::make_unique<t8_rbf_tps> ();
    }
    else {
      SC_ABORTF ("ERROR: RBF attribute missing or not correct. Unsupported RBF type.\n");
    }
  };

  /** Destructor. */
  ~t8_rbf () {};

  /**
   * Transfers the boundary node data to the internal data structure of the RBF.
   */
  void
  set_boundary_nodes (std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node>&& boundary_node_data)
  {
    boundary_nodes.clear ();
    boundary_nodes.reserve (boundary_node_data.size ());
    for (auto& [global_vertex_id, boundary_node] : boundary_node_data) {

      boundary_node.global_id = global_vertex_id;

      boundary_nodes.push_back (std::move (boundary_node));
    }
  }

  /** 
   * Search for the displacement of a specific boundary node with its associated global ID. 
   */
  t8_3D_vec
  get_boundary_displacement (const t8_gloidx_t global_id) const
  {
    for (const auto& node : boundary_nodes) {
      if (node.global_id == global_id) {
        return node.displacement;
      }
    }
    /* If the boundary node can not be found in the list of boundary nodes.*/
    t8_errorf ("ERROR: The boundary node %ld is missing in the boundary node list.\n", global_id);
    SC_ABORTF ("A boundary node is not recognized as one.");
  }

  /**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve ();

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

      const double psi = rbf_function->evaluate (distance, boundary_node.local_support_radius);
      /** Check if the basis function value (the influence factor) is not equal to zero with a numerical tolerance of 1e-12. */
      if (std::abs (psi) > 1e-12) {

        for (int coordinate = 0; coordinate < 3; ++coordinate) {
          /** Update the displacement for each coordinate with the weighted influence of the boundary node. */
          inner_node.displacement[coordinate] += boundary_node.weight[coordinate] * psi;
        }
      }
    }
  }

 private:
  /** The specialised solve function for either compactly supported or globally supported RBFs. 
   * \param[in] displacements The matrix of the boundary node displacements. Each row corresponds to a 
   *                          boundary node and the three columns correspond to the x, y, and z components of the displacement.
   * \param[in] num_boundary_nodes The number of boundary nodes.
   * \return The matrix of the calculated weights alpha.
  */
  Eigen::MatrixXd
  solve_compactly_supported_rbf (const Eigen::MatrixXd& displacements, const size_t num_boundary_nodes) const;
  Eigen::MatrixXd
  solve_globally_supported_rbf (const Eigen::MatrixXd& displacements, const size_t num_boundary_nodes) const;
  /** List of all registered support nodes. */
  std::vector<t8_rbf_boundary_node> boundary_nodes;
  /** Pointer to the radial basis function. */
  std::unique_ptr<t8_rbf_function> rbf_function;
};
