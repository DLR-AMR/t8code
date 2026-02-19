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

/** \file t8_cmesh_edge_conn_edge_to_tree.cxx
 *  This file implements the routines for the t8_cmesh_conn_edge_to_tree struct.
 */

#include <algorithm>
#include <memory>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_conn_edge_to_tree.hxx>

/* Builds edge_to_tree with existing tree_to_edge list. */
void
t8_cmesh_edge_conn_edge_to_tree::build_from_tte (const t8_cmesh_t cmesh, t8_cmesh_edge_conn_tree_to_edge& tte)
{
  T8_ASSERT (state == INITIALIZED);

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghosts = t8_cmesh_get_num_ghosts (cmesh);
  const t8_locidx_t num_local_trees_and_ghosts = num_local_trees + num_ghosts;

  for (t8_locidx_t itree = 0; itree < num_local_trees_and_ghosts; ++itree) {

    /* Get the global edge ids of this tree. */
    auto global_indices = tte.get_global_edges (cmesh, itree);

    /* Iterate over all local tree edges and add the global id to the list. */
    int iedge = 0;
    for (auto global_index : global_indices) {
      add_edge_to_tree (cmesh, global_index, itree, iedge++);
    }
  }

  /* Set state to committed. */
  state = COMMITTED;
}

/* Mark as ready for commit. Meaning that all
  * global edge ids have been added.
  * After commit, no edge ids can be added anymore. */
void
t8_cmesh_edge_conn_edge_to_tree::commit ([[maybe_unused]] const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  sort_list_by_tree_id ();
  state = COMMITTED;

  T8_ASSERT (contains_all_edges (cmesh));
}

void
t8_cmesh_edge_conn_edge_to_tree::add_edge_to_tree ([[maybe_unused]] const t8_cmesh_t cmesh,
                                                         t8_gloidx_t global_edge_id, t8_locidx_t ltreeid,
                                                         int tree_edge)
{
  T8_ASSERT (!is_committed ());
  T8_ASSERT (0 <= global_edge_id);
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
  const int num_tree_edges = t8_eclass_num_edges[tree_class];
  T8_ASSERT (0 <= tree_edge && tree_edge < num_tree_edges);
#endif
  if (state != INITIALIZED) {
    SC_ABORTF ("%s:%i: Trying to add edge to committed edge to tree structure.\n", __FILE__, __LINE__);
  }

  tree_edge_pair pair (ltreeid, tree_edge);
  tree_edge_list& list_of_globalid = edge_to_tree[global_edge_id];

  list_of_globalid.push_back (pair);
}

/** Compare function callback for two tree edge pairs. Implements the '<' operation for tree_edge_pair
 * (tree_id_A, edge_id_A) is consider smaller than
 * (tree_id_B, edge_id_B) if either
 *  tree_id_A < tree_id_B or
 *  tree_id_A == tree_id_B and edge_id_A < edge_id_B
 * \param [in] pair_a A pair of (local_tree_id_a, local_tree_edge_a)
 * \param [in] pair_b A pair of (local_tree_id_b, local_tree_edge_b)
 * \return True if  \a pair_a < \a pair_b meaning either \a local_tree_id_a < \a local_tree_id_b or
 *    (\a local_tree_id_a == \a local_tree_id_b and \a local_tree_edge_a < \a local_tree_edge_b)
 */
static inline bool
t8_cmesh_tree_edge_pair_compare (tree_edge_pair const& pair_a, tree_edge_pair const& pair_b)
{
  return pair_a.first == pair_b.first ?                              /* if tree_id_A == tree_id_B  */
           pair_a.second < pair_b.second                             /* then check edge_id_A < edge_id_B */
                                      : pair_a.first < pair_b.first; /* else check tree_id_A < tree_id_B */
}

void
t8_cmesh_edge_conn_edge_to_tree::sort_list_by_tree_id ()
{
  T8_ASSERT (!is_committed ());

  /* Iterate over each global edge */
  for (auto& [global_id, tree_edge_list] : edge_to_tree) {
    /* Check that the list contains at least one entry. */
    T8_ASSERT (tree_edge_list.size () > 0);
    /* Sort the list of local tree edges according to their
    * local tree id. */
    std::sort (tree_edge_list.begin (), tree_edge_list.end (), t8_cmesh_tree_edge_pair_compare);
  }
}

int
t8_cmesh_edge_conn_edge_to_tree::contains_all_edges (const t8_cmesh_t cmesh) const
{
  /* We need to check that each local tree/ghost and each edge
   * exists exactly once in the list.
   * We do so by setting up an indicator array storing the
   * number of edges for each tree and count down for each occurrence.
   * At the end the values must be zero. */

  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (is_committed ());

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);
  const t8_locidx_t num_trees_and_ghosts = num_local_trees + num_ghost_trees;

  std::vector<int> edge_counts (num_trees_and_ghosts);
  /* Fill each entry with the number of edges. */
  for (t8_locidx_t itree = 0; itree < num_trees_and_ghosts; ++itree) {
    /* Compute number of edges of this tree. */
    /* Get the trees class depending on whether it is a local tree or ghost. */
    const t8_eclass_t tree_class = itree < num_local_trees ? t8_cmesh_get_tree_class (cmesh, itree)
                                                           : t8_cmesh_get_ghost_class (cmesh, itree - num_local_trees);

    const int num_tree_edges = t8_eclass_num_edges[tree_class];

    /* Set the entry to the number of edges. */
    edge_counts[itree] = num_tree_edges;
  }

  /* Iterate over all entries in ett.
   * Each entry corresponds to a global edge id and
   * gives its list of tree indices and edges. */
  for (auto& [global_edge, tree_edge_list] : edge_to_tree) {
    /* Iterate over the list of tree indices and edges of this global edge. */
    for (auto& [tree_index, tree_edge] : tree_edge_list) {
      SC_CHECK_ABORT (0 <= tree_index && tree_index < num_trees_and_ghosts,
                      "Invalid tree id stored in edge to tree list.");
      const t8_eclass_t tree_class = tree_index < num_local_trees
                                       ? t8_cmesh_get_tree_class (cmesh, tree_index)
                                       : t8_cmesh_get_ghost_class (cmesh, tree_index - num_local_trees);

      const int num_tree_edges = t8_eclass_num_edges[tree_class];

      SC_CHECK_ABORT (0 <= tree_edge && tree_edge < num_tree_edges,
                      "Invalid edge id stored in edge to tree list.");

      /* remove this tree_edge from the edge_count */
      edge_counts[tree_index]--;
      /* Count must be >= 0 */
      T8_ASSERT (edge_counts[tree_index] >= 0);
    }
  }

  /* Now all entries must be set to 0 */
  for (int& entry : edge_counts) {
    if (entry != 0) {
      return 0;
    }
  }
  return 1;
}
