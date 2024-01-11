/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_vertex_conn_tree_to_vertex_c.hxx
 * Class to save data structure for tree_to_vertex_lists of the cmesh.
 * When the cmesh stores global vertex numbers, we require a lookup that 
 * matches a tree and its local vertex to a global vertex id.
 * This lookup is encoded in the t8_cmesh_vertex_conn_tree_to_vertex_c struct.
 */

/* TODO:
 *  It is probably best to set all global ids of a single tree as one attribute.
 * That way we can store it as a single arrays of id's.
 * We will probably most often want to access all ids of one tree, so a function returning an array of ids
 * is a good idea anyway and if we already store them as such then we do not need to do any data movement
 * when accessing.
 * 
 * On the downside we will only have a "set all ids of a tree" function and no "set this single id for this tree and this vertex" function.
 */

#ifndef T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX
#define T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX

#include <t8_cmesh.h>

typedef struct t8_cmesh_vertex_conn_tree_to_vertex
{
 public:
  /* Setter functions */
  /** Set the global vertex id of a local vertex of a local tree.
   * \param[in] cmesh The considered cmesh
   * \param[in] local_tree A local tree id of \a cmesh
   * \param[in] local_tree_vertex A vertex of \a local_tree. 0 <= \a local_tree_vertex < NUM_VERTICES(\a local_tree)
   * \param[in] global_vertex_id The id of the global vertex.
   * 
   * \note Cmesh must not be committed.
  */
  void
  set_global_vertex_id_of_tree_vertex (const t8_cmesh_t, const t8_locidx_t local_tree, const int local_tree_vertex,
                                       const int global_vertex_id);

  /** Set all global vertex ids of a local tree. 
   * \param[in] cmesh The considered cmesh
   * \param[in] local_tree A local tree id of \a cmesh
   * \param[in] global_vertex_id The ids of the global vertices in order of \a local_tree's vertices.
   * \param[in] num_vertices Must match the number of vertices of \a local_tree
   * 
   * \note Cmesh must not be committed.
  */
  void
  set_global_vertex_ids_of_tree_vertices (const t8_cmesh_t, const t8_locidx_t local_tree, const int *global_tree_vertex,
                                          const int num_vertices);

  t8_gloidx_t
  get_global_vertex (const t8_cmesh_t, const t8_locidx_t local_tree, const int local_tree_vertex);

  const t8_gloidx_t *
  get_global_vertices (const t8_cmesh_t, const t8_locidx_t local_tree, const int num_vertices);

 private:
} t8_cmesh_vertex_conn_tree_to_vertex_c;

#endif /* !T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX */
