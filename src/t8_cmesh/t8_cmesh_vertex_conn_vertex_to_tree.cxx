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

#include <stdexcept>
#include <algorithm>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_conn_vertex_to_tree.hxx>

/** \file t8_cmesh_conn_vertex_to_tree.cxx
 *  This file implements the routines for the t8_cmesh_conn_vertex_to_tree_c struct.
 */

t8_cmesh_tree_vertex_list&
t8_cmesh_vertex_conn_vertex_to_tree_c::get_tree_list_of_vertex (t8_gloidx_t global_vertex_id)
{
  T8_ASSERT (is_committed ());
  T8_ASSERT (0 <= global_vertex_id);

  /* Use at() to look for the vertex entry.
   * If the entry does not exist an exception of
   * type std::out_of_range is thrown. */

  try {
    return vertex_to_tree.at (global_vertex_id);
  } catch (const std::out_of_range& e) {
    t8_errorf ("ERROR: Could not find vertex %li for cmesh.\n", global_vertex_id);
    SC_ABORTF ("Caught exception 'out of range': %s\n", e.what ());
  }
}

int
t8_cmesh_vertex_conn_vertex_to_tree_c::is_committed ()
{
  return state == COMMITTED;
}

/* Mark as ready for commit. Meaning that all 
  * global vertex ids have been added.
  * After commit, no vertex ids can be added anymore. */
void
t8_cmesh_vertex_conn_vertex_to_tree_c::commit (t8_cmesh_t cmesh)
{
  /* TODO: In debugging mode, check whether all local trees have a global id. */
  /* Sort the entries of the global ids. */
  sort_list_by_tree_id ();
  state = COMMITTED;
}

void
t8_cmesh_vertex_conn_vertex_to_tree_c::add_vertex_to_tree (t8_cmesh_t cmesh, t8_gloidx_t global_vertex_id,
                                                           t8_locidx_t ltreeid, int tree_vertex)
{
  T8_ASSERT (!is_committed ());
  T8_ASSERT (0 <= global_vertex_id);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid) || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

#if T8_ENABLE_DEBUG
  t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  int num_tree_vertices = t8_eclass_num_vertices[tree_class];
  T8_ASSERT (0 <= tree_vertex && tree_vertex < num_tree_vertices);
#endif
  if (state != INITIALIZED) {
    SC_ABORTF ("Trying to add vertex to committed vertex to tree structure.\n");
  }

  t8_cmesh_tree_vertex_pair pair (ltreeid, tree_vertex);
  t8_cmesh_tree_vertex_list& list_of_globalid = vertex_to_tree[global_vertex_id];

  list_of_globalid.push_back (pair);
}

/* Compare two tree vertex pairs. 
 * (tree_id_A, vertex_id_A) is consider smaller than
 * (tree_id_B, vertex_id_B) if either
 *  tree_id_A < tree_id_B or
 *  tree_id_A == tree_id_B and vertex_id_A < vertex_id_B */
static int
t8_cmesh_tree_vertex_pair_compare (t8_cmesh_tree_vertex_pair const& pair_a, t8_cmesh_tree_vertex_pair const& pair_b)
{
  return pair_a.first == pair_b.first ?                              /* if tree_id_A == tree_id_B  */
           pair_a.second < pair_b.second                             /* then check vertex_id_A < vertex_id_B */
                                      : pair_a.first < pair_b.first; /* else check tree_id_A < tree_id_B */
}

void
t8_cmesh_vertex_conn_vertex_to_tree_c::sort_list_by_tree_id ()
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