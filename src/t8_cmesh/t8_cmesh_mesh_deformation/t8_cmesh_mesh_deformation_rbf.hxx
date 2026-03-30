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

/** \file t8_cmesh_mesh_deformation_rbf.hxx
 * Implementation of CAD-based mesh deformation.
 */
#pragma once

#include <vector>
#include <cmath>
#include <t8_types/t8_vec.hxx>

/** 
 * The available RBF function types. 
 */
typedef enum { T8_RBF_CP_C2 = 0, T8_RBF_TPS, T8_RBF_FUNCTION_COUNT } t8_rbf_function_type;

/**
 * This struct will be a single boundary node from the mesh deformation. 
 * The new position of this node can be directly calculated from the new incoming CAD geometry given. 
 */
struct t8_rbf_boundary_node
{
  /* The old position of the boundary node. */
  t8_3D_vec position;
  /* The displacement known and calculated from the different coordinates of the CAD geometries. */
  t8_3D_vec displacement;
  /* The calculated RBF coefficient alpha. */
  t8_3D_vec weight;
};

/**  
 * Struct for mesh deformation using Radial Basis Functions (RBF).
 * It handles the interpolation of the boundary displacements to the inner nodes. 
 */
struct t8_cmesh_mesh_deformation_rbf
{
 public:
  /** Constructor. 
   * \param[in] support_radius  The radius of the compact local support for the Wendland CP C2 function.
   *                            If the TPS RBF is used, we will not use the support_radius due to its global support.
   * \param[in] rbf_type        The RBF type to be used.
  */
  t8_cmesh_mesh_deformation_rbf (double support_radius, t8_rbf_function_type rbf_type)
    : radius (support_radius), selected_type (rbf_type) {};

  /** Destructor. */
  ~t8_cmesh_mesh_deformation_rbf () {};

  /**
   * Add a new support node to the RBF system based on the boundary node  (they are the same if we do not use a greedy algorithm to get less support nodes).
   * These nodes are then used in the linear system A * alpha = d.
   * The coefficient alpha is the needed weight of a support node which ensures that the RBF interpolation 
   * exactly matches the CAD displacement because it corrects the spatial overlap of nearby basis functions.
   * \param[in] position      The coordinates of the node before the displement.
   * \param[in] displacement  The known displacement from the input geometries.
   */
  void
  add_node (const double position[3], const double displacement[3]);

  /**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
  void
  solve ();

  /**
   * Interpolates the inner node.
   * \param[in] inner_node
   * \param[out] inner_node_displacement
   */
  void
  interpolate (const double inner_node[3], const double inner_node_displacement[3]);

 private:
  /** Wendland function. psi(x) = (1 - x)_+^4 * (4x + 1) for (1-x) > 0. 
   * \param[in] distance    The euclidean distance between two points.
   * return     The computed weight (when using solve()) or influence (when using interpolate()).
  */
  const double
  wendland_cp_c2 (const double distance);

  /** Thin Plate Spline: psi(x) = x^2 * log(x). 
   * \param[in] distance  The euclidean distance between two points.
   * return     The computed weight (when using solve()) or influence (when using interpolate()).
  */
  const double
  thin_plate_spline (const double distance);
  /** The chosen radius for the local support of the RBF CP C^2. */
  double radius;
  /** List of all registered support nodes. */
  std::vector<t8_rbf_boundary_node> boundary_nodes;
  /** Selected Radial Basis Function. */
  t8_rbf_function_type selected_type;
};
