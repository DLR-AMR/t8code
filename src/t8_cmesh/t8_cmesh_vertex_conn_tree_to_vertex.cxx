/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_conn_tree_to_vertex.hxx>

/** \file t8_cmesh_conn_tree_to_vertex.cxx
 *  This file implements the routines for the t8_cmesh_conn_tree_to_vertex_c struct.
 */

/* Set all global vertex ids of a local tree. 
* \param[in] cmesh The considered cmesh
* \param[in] local_tree A local tree id of \a cmesh
* \param[in] global_vertex_id The ids of the global vertices in order of \a local_tree's vertices.
* \param[in] num_vertices Must match the number of vertices of \a local_tree
* 
* \note Cmesh must not be committed.
*/
void
t8_cmesh_vertex_conn_tree_to_vertex::set_global_vertex_ids_of_tree_vertices (const t8_cmesh_t cmesh,
                                                                             const t8_gloidx_t global_tree,
                                                                             const t8_gloidx_t *global_tree_vertices,
                                                                             const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (num_vertices >= 0);
  T8_ASSERT (global_tree_vertices != NULL);

  /* TODO: we currently do not check whether the num_vertices argument
   *       matches the number of vertices of the tree.
   *       We cannot do it here, since this function call happens before commit,
   *       thus we might not even know the eclass of the tree.
   *       Maybe it is possible to check this during t8_cmesh_add_attributes?
   */

  /* We copy the data directly, hence set data_persiss to 0 */
  const int data_persists = 0;
  t8_debugf ("Setting %i global vertices for global tree %li.\n", num_vertices, global_tree);
  t8_cmesh_set_attribute_gloidx_array (cmesh, global_tree, t8_get_package_id (), T8_CMESH_GLOBAL_VERTICES_ATTRIBUTE_KEY,
                                       global_tree_vertices, num_vertices, data_persists);
}

/* TODO: What if the attribute is not set? error handling */
const t8_gloidx_t *
t8_cmesh_vertex_conn_tree_to_vertex::get_global_vertices (const t8_cmesh_t cmesh, const t8_locidx_t local_tree,
                                                          const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

#if T8_ENABLE_DEBUG
  /* Verify that num_vertices matches the number of tree vertices */
  const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, local_tree);
  const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

  T8_ASSERT (num_vertices == num_tree_vertices);
#endif

  t8_debugf ("Getting %i global vertices for local tree %i.\n", num_vertices, local_tree);
  const t8_gloidx_t *global_vertices = t8_cmesh_get_attribute_gloidx_array (
    cmesh, t8_get_package_id (), T8_CMESH_GLOBAL_VERTICES_ATTRIBUTE_KEY, local_tree, num_vertices);
  T8_ASSERT (global_vertices != NULL);
  return global_vertices;
}

/* TODO: What if the attribute is not set? error handling */
t8_gloidx_t
t8_cmesh_vertex_conn_tree_to_vertex::get_global_vertex (const t8_cmesh_t cmesh, const t8_locidx_t local_tree,
                                                        const int local_tree_vertex, const int num_tree_vertices)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  /* Verify that local_tree_vertex is in fact a local vertex of the tree */
  /* Note: We only perform this check in debugging mode.
   *       In non-debugging mode, using a vertex index beyond the trees index allows
   *       for a potential attacker to gain access to memory possibly not owned by the caller.
   *       We do not check in non-debugging mode for (obvious) performance reasons. */
  T8_ASSERT (0 <= local_tree_vertex);
  T8_ASSERT (local_tree_vertex < num_tree_vertices);

  return get_global_vertices (cmesh, local_tree, num_tree_vertices)[local_tree_vertex];
}