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

/** \file t8_cmesh_conn_tree_to_vertex.cxx
 *  This file implements the routines for the t8_cmesh_conn_tree_to_vertex struct.
 */

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.hxx>

/* constructor from a given vertex to tree list. */
t8_cmesh_vertex_conn_tree_to_vertex::t8_cmesh_vertex_conn_tree_to_vertex (
  const t8_cmesh_t cmesh_from, const t8_cmesh_t cmesh, const t8_cmesh_vertex_conn_vertex_to_tree &vtt)
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

  /* Stores an array of local indices and mapped global vertices for each tree.
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

