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

/** \file t8_cmesh_vertex_connectivity.h
 * C interface for the cmesh vertex connectivity. See \ref t8_cmesh_vertex_connectivity.hxx
 * for further information.
 */

#ifndef T8_CMESH_VERTEX_CONNECTIVITY_H
#define T8_CMESH_VERTEX_CONNECTIVITY_H

T8_EXTERN_C_BEGIN ();

typedef struct t8_cmesh_vertex_connectivity *t8_cmesh_vertex_connectivity_c;

/** Set all global vertex ids of a local tree.
  * \param [in] cmesh The considered cmesh
  * \param [in] global_tree A global tree id of \a cmesh
  * \param [in] global_tree_vertices The ids of the global vertices in order of \a local_tree's vertices.
  * \param [in] num_vertices Must match the number of vertices of \a local_tree
  *
  * \note \a cmesh must not be committed.
  */
void
t8_cmesh_set_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                      const t8_gloidx_t *global_tree_vertices, const int num_vertices);

/** Return the total number of global vertices of a cmesh (across all processes).
 * \return The total number of global vertices of \a cmesh.
 * \note \a cmesh must be committed.
*/
t8_gloidx_t
t8_cmesh_get_num_global_vertices (const t8_cmesh_t cmesh);

/** Return the number of process local global vertices of a cmesh.
 * \return The number of process local global vertices of \a cmesh.
 * \note \a cmesh must be committed.
*/
t8_locidx_t
t8_cmesh_get_num_local_vertices (const t8_cmesh_t cmesh);

/** Get the global vertex indices of a tree in its local vertex order.
 * \param [in] cmesh A committed cmesh.
 * \param [in] local_tree A local tree in \a cmesh.
 * \param [in] num_vertices The number of vertices of \a local_tree
 * \return The global vertices of \a local_tree
 */
const t8_gloidx_t *
t8_cmesh_get_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int num_vertices);

/** Get a single global vertex index of a local tree's local vertex.
* \param [in] cmesh A committed cmesh.
* \param [in] local_tree A local tree in \a cmesh.
* \param [in] local_tree_vertex A local vertex of \a local_tree
* \param [in] num_vertices The number of vertices of \a local_tree
* \return The global vertex matching \a local_tree_vertex of \a local_tree.
*/
t8_gloidx_t
t8_cmesh_get_global_vertex_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_vertex,
                                    const int num_vertices);
/** Get the number of global trees a global vertex is connected to.
 * \param [in] global_vertex_id The global id of a vertex in the cmesh.
 * \note if a tree is contained multiple times it is counted as multiple entries.
 * Example: For a quad where all 4 vertices map to a single global vertex this function will return 4.
 */
int
t8_cmesh_get_num_trees_at_vertex (const t8_cmesh_t cmesh, t8_gloidx_t global_vertex);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_VERTEX_CONNECTIVITY_H */
