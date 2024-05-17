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

/* constructor from a given vertex to tree list. */
t8_cmesh_vertex_conn_tree_to_vertex::t8_cmesh_vertex_conn_tree_to_vertex (
  const t8_cmesh_t cmesh_from, t8_cmesh_t cmesh, const t8_cmesh_vertex_conn_vertex_to_tree &vtt)
  : t8_cmesh_vertex_conn_tree_to_vertex ()
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh_from));
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (vtt.is_committed ());

  /* Since we need to collect all global ids of a tree before we can add the tree,
   * we need a temporary buffer structure.
   * Note that this is very memory intensive, doubling the required memory. */

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh_from);
  const t8_locidx_t num_ghosts = t8_cmesh_get_num_ghosts (cmesh_from);
  const t8_locidx_t num_local_trees_and_ghosts = num_local_trees + num_ghosts;

  /* Stores for each tree an array of local index and mapped global vertex.
   * We need to store the local indices as well, since they will not be
   * sorted. Hence we need to sort later. */
  using global_local_index_pair = std::pair<t8_gloidx_t, t8_locidx_t>;
  std::vector<std::vector<global_local_index_pair>> global_ids_per_tree (num_local_trees_and_ghosts);

  /* Iterate over all entries of the vertex to tree list. */
  for (auto &[global_vertex, vertex_list] : vtt) {
    /* Iterate over all entries of the vertex list,
     * thus, each entry is a local tree id and vertex index. */
    for (auto &[tree_index, vertex_index] : vertex_list) {
      global_local_index_pair new_pair (global_vertex, vertex_index);
      global_ids_per_tree[tree_index].push_back (new_pair);
    }
  }
  /* Now sort the list and add the global indices of each tree. */
  for (t8_locidx_t itree = 0; itree < num_local_trees_and_ghosts; ++itree) {
    auto &tree_indices = global_ids_per_tree[itree];
    /* Sort according to local index. */
    std::sort (tree_indices.begin (), tree_indices.end (),
               [] (global_local_index_pair &a, global_local_index_pair &b) { return a.second < b.second; });

    /* Compute number of vertices of this tree. */
    /* Get the trees class depending on whether it is a local tree or ghost. */
    const t8_eclass_t tree_class = itree < num_local_trees
                                     ? t8_cmesh_get_tree_class (cmesh_from, itree)
                                     : t8_cmesh_get_ghost_class (cmesh_from, itree - num_local_trees);

    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Size of list must match number of tree vertices. */
    SC_CHECK_ABORTF (
      tree_indices.size () == (size_t) num_tree_vertices,
      "ERROR. Number of mapped local tree vertices for tree %i does not equal number of tree vertices. %zu != %i\n",
      itree, tree_indices.size (), num_tree_vertices);
    /* We now build the array of global indices of this tree, sorted 
     * by the trees local indices. */
    t8_gloidx_t *global_tree_indices = T8_ALLOC (t8_gloidx_t, num_tree_vertices);
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      global_tree_indices[ivertex] = tree_indices[ivertex].first;
    }
    const t8_gloidx_t global_tree_id = t8_cmesh_get_global_id (cmesh_from, itree);
    set_global_vertex_ids_of_tree_vertices (cmesh, global_tree_id, global_tree_indices, num_tree_vertices);
    T8_FREE (global_tree_indices);
  }
}

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
                                                          const int num_vertices) const
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
                                                        const int local_tree_vertex, const int num_tree_vertices) const
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