
/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#ifndef T8_CMESH_DEFORMATION_HXX
#define T8_CMESH_DEFORMATION_HXX

#include <memory>
#include <unordered_map>
#include <array>

#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cad/t8_cad.hxx>
#include <t8_types/t8_vec.hxx>

/** Struct for mesh deformation. */
struct t8_cmesh_mesh_deformation
{
 public:
  /** Constructor */
  t8_cmesh_mesh_deformation (): associated_cmesh (nullptr), updated_geometry (nullptr) {};

  /** Destructor */
  ~t8_cmesh_mesh_deformation () {};

 private:
  /** A pointer to the cmesh for attribute retrieval */
  t8_cmesh_t associated_cmesh;
  /** A shared pointer to the updated geometry which comes from a new cad file */
  std::shared_ptr<t8_cad> updated_geometry;
};

/** 
 * Computes the displacements of the surface vertices.
 * 
 * \param [in] cmesh The committed mesh
 * \param [in] cad Pointer to the CAD data
 * \return Map from global vertex ID to 3D displacement vector
 */
std::unordered_map<t8_gloidx_t, t8_3D_vec>
calculate_displacement_surface_vertices (t8_cmesh_t cmesh, const t8_cad *cad);

/**
 * Apply vertex displacements to a committed cmesh.
 *
 * Iterates over the provided map of global vertex IDs to 3D displacement vectors,
 * updating the coordinates in each tree where the vertex appears.
 *
 * @param cmesh The committed coarse mesh structure.
 * @param displacements Map from global vertex ID to 3D displacement vector [dx, dy, dz].
 */
void
apply_vertex_displacements (t8_cmesh_t cmesh, const std::unordered_map<t8_gloidx_t, t8_3D_vec> &displacements,
                            std::shared_ptr<t8_cad> cad);
#endif /* !T8_CMESH_DEFORMATION_HXX */
