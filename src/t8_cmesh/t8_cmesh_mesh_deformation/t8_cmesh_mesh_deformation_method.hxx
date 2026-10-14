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

/** \file t8_cmesh_mesh_deformation_method.hxx
 *  This file implements the base class for the mesh deformation methods.
 */

#pragma once

#include <t8.h>
#include <t8_types/t8_vec.hxx>
#include <unordered_map>

/** A struct representing a node in the mesh for the deformation, containing its global ID, position, displacement and the mean distance to its neighbors. */
struct t8_mesh_deformation_node
{
  /**The global ID of the node. */
  t8_gloidx_t global_id;
  /** The old position of the node. */
  t8_3D_vec position;
  /** The displacement. 
   * For boundary nodes this is known and calculated from the different coordinates of the CAD geometries. 
   * For inner nodes this is the result of the interpolation. */
  t8_3D_vec displacement;
  /** Mean distance to neighbouring vertices (boundary nodes only). */
  double mean_neighbor_distance = 0.0;
};

/**
 * The base class for all mesh deformation methods.
 * This class provides a general template for all methods.
 * It is a pure virtual class and has to be inherited by a concrete
 * mesh deformation method.
 */
struct t8_mesh_deformation_method
{
  virtual ~t8_mesh_deformation_method () = default;

  /**
   * Transfers the boundary node data to the internal data structure of the mesh defromation.
   * \param [in] boundary_node_data Map of global vertex IDs to boundary nodes.
   */
  virtual void
  set_boundary_nodes (const std::unordered_map<t8_gloidx_t, t8_mesh_deformation_node>& boundary_node_data)
    = 0;

  /**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory for RBF to later be able to interpolate the inner nodes.
   */
  virtual void
  solve ()
  {
  }

  /**
   * Interpolates the inner node.
   * \param [in, out] inner_node   The inner node which will be interpolated.
   */
  virtual void
  interpolate (t8_mesh_deformation_node& inner_node) const
    = 0;
};
