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

/** \file t8_cmesh_edge_connectivity.h
 * C interface for the cmesh edge connectivity. See \ref t8_cmesh_edge_connectivity.hxx
 * for further information.
 */

#ifndef T8_CMESH_EDGE_CONNECTIVITY_H
#define T8_CMESH_EDGE_CONNECTIVITY_H

T8_EXTERN_C_BEGIN ();

/**
 * t8_cmesh_edge_connectivity_c
 *
 * Opaque pointer to the cmesh edge connectivity structure.
 */
typedef struct t8_cmesh_edge_connectivity *t8_cmesh_edge_connectivity_c;

/** Set all global edge ids of a local tree.
  * \param [in] cmesh The considered cmesh
  * \param [in] global_tree A global tree id of \a cmesh
  * \param [in] global_tree_edges The ids of the global edges in order of \a local_tree's edges.
  * \param [in] num_edges Must match the number of edges of \a local_tree
  *
  * \note \a cmesh must not be committed.
  */
void
t8_cmesh_set_global_edges_of_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                      const t8_gloidx_t *global_tree_edges, const int num_edges);

/** Return the total number of global edges of a cmesh (across all processes).
 * \return The total number of global edges of \a cmesh.
 * \note \a cmesh must be committed.
*/
t8_gloidx_t
t8_cmesh_get_num_global_edges (const t8_cmesh_t cmesh);

/** Return the number of process local global edges of a cmesh.
 * \return The number of process local global edges of \a cmesh.
 * \note \a cmesh must be committed.
*/
t8_locidx_t
t8_cmesh_get_num_local_edges (const t8_cmesh_t cmesh);

/** Get the global edge indices of a tree in its local edge order.
 * \param [in] cmesh A committed cmesh.
 * \param [in] local_tree A local tree in \a cmesh.
 * \param [out] num_edges The number of edges of \a local_tree, if not null
 * \return The global edges of \a local_tree
 */
const t8_gloidx_t *
t8_cmesh_get_global_edges_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, int *num_edges);

/** Get a single global edge index of a local tree's local edge.
* \param [in] cmesh A committed cmesh.
* \param [in] local_tree A local tree in \a cmesh.
* \param [in] local_tree_edge A local edge of \a local_tree
* \return The global edge matching \a local_tree_edge of \a local_tree.
*/
t8_gloidx_t
t8_cmesh_get_global_edge_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_edge);

/** Get the number of global trees a global edge is connected to.
 * \param [in] cmesh A committed cmesh.
 * \param [in] global_edge The global id of a edge in the cmesh.
 * \note if a tree is contained multiple times it is counted as multiple entries.
 * Example: For a quad where all 4 edges map to a single global edge this function will return 4.
 */
int
t8_cmesh_get_num_trees_at_edge (const t8_cmesh_t cmesh, t8_gloidx_t global_edge);

/** Check whether a given cmesh uses the edge connectivity feature.
 * \param [in] cmesh A committed cmesh.
 * \return Nonzero if the \a cmesh uses the edge connectivity feature, zero if not.
 */
int
t8_cmesh_uses_edge_connectivity (const t8_cmesh_t cmesh);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_EDGE_CONNECTIVITY_H */
