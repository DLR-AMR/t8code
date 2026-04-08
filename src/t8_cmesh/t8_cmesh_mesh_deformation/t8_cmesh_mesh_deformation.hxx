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

/** \file t8_cmesh_mesh_deformation.hxx
 * Implementation of CAD-based mesh deformation.
 */

#pragma once

#include <memory>
#include <unordered_map>
#include <array>

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cad/t8_cad_handle.hxx>
#include <t8_types/t8_vec.hxx>

#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_rbf.hxx>

/** Struct for mesh deformation. */
struct t8_cmesh_mesh_deformation
{
 public:
  /** Constructor */
  t8_cmesh_mesh_deformation (t8_cmesh_t cmesh): associated_cmesh (cmesh), updated_geometry (nullptr) {};

  /** Destructor */
  ~t8_cmesh_mesh_deformation () {};

  /** 
 * Computes the displacements of the surface vertices.
 * 
 * \param [in]  cad A pointer to the CAD-based geometry object.
 * \return Map from global vertex ID to RBF boundary node which contains the displacement and can than be used to calculate the weight of the boundary node.
 */
  std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node>
  calculate_displacement_surface_vertices (const t8_cad_handle *cad);

  /**
 * Apply vertex displacements to a committed cmesh.
 *
 * Iterates over the provided map of global vertex IDs to 3D displacement vectors,
 * updating the coordinates in each tree where the vertex appears.
 *
 * \param [in] boundary_node_data Map from global vertex ID to RBF boundary node.
 * \param [in] cad The shared pointer to the CAD geometry to update.
 */
  void
  apply_vertex_displacements (const std::unordered_map<t8_gloidx_t, t8_rbf_boundary_node> &boundary_node_data,
                              std::shared_ptr<t8_cad_handle> cad);

 private:
  /** A pointer to the cmesh for attribute retrieval */
  t8_cmesh_t associated_cmesh;
  /** A shared pointer to the updated geometry which comes from a new cad file */
  std::shared_ptr<t8_cad_handle> updated_geometry;
};
