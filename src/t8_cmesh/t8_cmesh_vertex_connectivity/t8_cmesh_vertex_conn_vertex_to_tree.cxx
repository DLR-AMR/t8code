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

/** \file t8_cmesh_conn_vertex_to_tree.cxx
 *  This file implements the routines for the t8_cmesh_conn_vertex_to_tree struct.
 */

#include <algorithm>
#include <memory>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.hxx>

/* Builds vertex_to_tree with existing tree_to_vertex list. */
void
t8_cmesh_vertex_conn_vertex_to_tree::build_from_ttv (const t8_cmesh_t cmesh, t8_cmesh_vertex_conn_tree_to_vertex& ttv)
{
  T8_ASSERT (state == INITIALIZED);

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghosts = t8_cmesh_get_num_ghosts (cmesh);
  const t8_locidx_t num_local_trees_and_ghosts = num_local_trees + num_ghosts;

  for (t8_locidx_t itree = 0; itree < num_local_trees_and_ghosts; ++itree) {
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Get the global vertex ids of this tree. */
    const t8_gloidx_t* global_indices = ttv.get_global_vertices (cmesh, itree, num_tree_vertices);

    /* Iterate over all local tree vertices and add the global id to the list. */
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      add_vertex_to_tree (cmesh, global_indices[ivertex], itree, ivertex);
    }
  }

  /* Set state to committed. */
  state = COMMITTED;
}

/* Mark as ready for commit. Meaning that all 
  * global vertex ids have been added.
  * After commit, no vertex ids can be added anymore. */
void
t8_cmesh_vertex_conn_vertex_to_tree::commit ([[maybe_unused]] const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  sort_list_by_tree_id ();
  state = COMMITTED;

  T8_ASSERT (contains_all_vertices (cmesh));
}

void
t8_cmesh_vertex_conn_vertex_to_tree::add_vertex_to_tree ([[maybe_unused]] const t8_cmesh_t cmesh,
                                                         t8_gloidx_t global_vertex_id, t8_locidx_t ltreeid,
                                                         int tree_vertex)
{
  T8_ASSERT (!is_committed ());
  T8_ASSERT (0 <= global_vertex_id);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid) || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

#if T8_ENABLE_DEBUG
  const bool is_local_tree = t8_cmesh_treeid_is_local_tree (cmesh, ltreeid);
  t8_eclass_t tree_class;
  if (is_local_tree) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  }
  else {
    tree_class = t8_cmesh_get_ghost_class (cmesh, ltreeid - t8_cmesh_get_num_local_trees (cmesh));
  }
  const int num_tree_vertices = t8_eclass_num_vertices[tree_class];
  T8_ASSERT (0 <= tree_vertex && tree_vertex < num_tree_vertices);
#endif
  if (state != INITIALIZED) {
    SC_ABORTF ("%s:%i: Trying to add vertex to committed vertex to tree structure.\n", __FILE__, __LINE__);
  }

  tree_vertex_pair pair (ltreeid, tree_vertex);
  tree_vertex_list& list_of_globalid = vertex_to_tree[global_vertex_id];

  list_of_globalid.push_back (pair);
}

/** Compare function callback for two tree vertex pairs. Implements the '<' operation for tree_vertex_pair
 * (tree_id_A, vertex_id_A) is consider smaller than
 * (tree_id_B, vertex_id_B) if either
 *  tree_id_A < tree_id_B or
 *  tree_id_A == tree_id_B and vertex_id_A < vertex_id_B
 * \param [in] pair_a A pair of (local_tree_id_a, local_tree_vertex_a)
 * \param [in] pair_b A pair of (local_tree_id_b, local_tree_vertex_b)
 * \return True if  \a pair_a < \a pair_b meaning either \a local_tree_id_a < \a local_tree_id_b or
 *    (\a local_tree_id_a == \a local_tree_id_b and \a local_tree_vertex_a < \a local_tree_vertex_b)
 */
static inline bool
t8_cmesh_tree_vertex_pair_compare (tree_vertex_pair const& pair_a, tree_vertex_pair const& pair_b)
{
  return pair_a.first == pair_b.first ?                              /* if tree_id_A == tree_id_B  */
           pair_a.second < pair_b.second                             /* then check vertex_id_A < vertex_id_B */
                                      : pair_a.first < pair_b.first; /* else check tree_id_A < tree_id_B */
}

void
t8_cmesh_vertex_conn_vertex_to_tree::sort_list_by_tree_id ()
{
  T8_ASSERT (!is_committed ());

  /* Iterate over each global vertex */
  for (auto& [global_id, tree_vertex_list] : vertex_to_tree) {
    /* Check that the list contains at least one entry. */
    T8_ASSERT (tree_vertex_list.size () > 0);
    /* Sort the list of local tree vertices according to their 
    * local tree id. */
    std::sort (tree_vertex_list.begin (), tree_vertex_list.end (), t8_cmesh_tree_vertex_pair_compare);
  }
}

int
t8_cmesh_vertex_conn_vertex_to_tree::contains_all_vertices (const t8_cmesh_t cmesh) const
{
  /* We need to check that each local tree/ghost and each vertex 
   * exists exactly once in the list. 
   * We do so by setting up an indicator array storing the
   * number of vertices for each tree and count down for each occurrence.
   * At the end the values must be zero. */

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (is_committed ());

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);
  const t8_locidx_t num_trees_and_ghosts = num_local_trees + num_ghost_trees;

  std::vector<int> vertex_counts (num_trees_and_ghosts);
  /* Fill each entry with the number of vertices. */
  for (t8_locidx_t itree = 0; itree < num_trees_and_ghosts; ++itree) {
    /* Compute number of vertices of this tree. */
    /* Get the trees class depending on whether it is a local tree or ghost. */
    const t8_eclass_t tree_class = itree < num_local_trees ? t8_cmesh_get_tree_class (cmesh, itree)
                                                           : t8_cmesh_get_ghost_class (cmesh, itree - num_local_trees);

    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Set the entry to the number of vertices. */
    vertex_counts[itree] = num_tree_vertices;
  }

  /* Iterate over all entries in vtt.
   * Each entry corresponds to a global vertex id and
   * gives its list of tree indices and vertices. */
  for (auto& [global_vertex, tree_vertex_list] : vertex_to_tree) {
    /* Iterate over the list of tree indices and vertices of this global vertex. */
    for (auto& [tree_index, tree_vertex] : tree_vertex_list) {
      SC_CHECK_ABORT (0 <= tree_index && tree_index < num_trees_and_ghosts,
                      "Invalid tree id stored in vertex to tree list.");
      const t8_eclass_t tree_class = tree_index < num_local_trees
                                       ? t8_cmesh_get_tree_class (cmesh, tree_index)
                                       : t8_cmesh_get_ghost_class (cmesh, tree_index - num_local_trees);

      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

      SC_CHECK_ABORT (0 <= tree_vertex && tree_vertex < num_tree_vertices,
                      "Invalid vertex id stored in vertex to tree list.");

      /* remove this tree_vertex from the vertex_count */
      vertex_counts[tree_index]--;
      /* Count must be >= 0 */
      T8_ASSERT (vertex_counts[tree_index] >= 0);
    }
  }

  /* Now all entries must be set to 0 */
  for (int& entry : vertex_counts) {
    if (entry != 0) {
      return 0;
    }
  }
  return 1;
}
