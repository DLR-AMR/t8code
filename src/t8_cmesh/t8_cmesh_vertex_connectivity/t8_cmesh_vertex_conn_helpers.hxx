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

/** \file t8_cmesh_vertex_conn_helpers.hxx
 * Helper routines for or using the vertex connectivity.
 */

#ifndef T8_CMESH_VERTEX_CONN_HELPERS_HXX
#define T8_CMESH_VERTEX_CONN_HELPERS_HXX

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

#include <array>
#include <vector>

/** Derive global vertex ids for every tree currently in \a cmesh's stash
 * from the trees corner coordinates, and set them via
 * \ref t8_cmesh_set_global_vertices_of_tree, once per tree. Tree corners
 * that coincide geometrically (within a tolerance derived from the mesh's
 * bounding box) get the same global vertex id.
 *
 * \param [in,out] cmesh  An initialized, uncommitted cmesh whose stash still
 *                        holds tree classes and vertex coordinates.
 *
 * \note Must be called before the stash is converted into the committed
 *       tree structure. Is therefore called during \ref t8_cmesh_commit.
 */
void
t8_cmesh_vertex_conn_set_vertices_by_coordinates (const t8_cmesh_t cmesh);

/** Cluster the corners of every tree currently in \a stash by coincident coordinates (within a
 * tolerance derived from the mesh's bounding box) and assign each cluster a consecutive global vertex
 * id, starting at 0.
 * \param [in] stash      The stash.
 * \param [in] eclasses   The eclasses of the stash as obtained by \ref t8_stash_extract_eclasses.
 * \return Global vertex ids, one fixed-size corner array per tree.
 */
std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>>
t8_cmesh_vertex_conn_helpers_global_ids_by_coordinates (const t8_stash_t stash,
                                                        const std::vector<t8_eclass_t> &eclasses);

#endif /* !T8_CMESH_VERTEX_CONN_HELPERS_HXX */
