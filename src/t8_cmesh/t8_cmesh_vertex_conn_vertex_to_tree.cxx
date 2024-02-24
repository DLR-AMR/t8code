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
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_conn_vertex_to_tree.hxx>

/** \file t8_cmesh_conn_vertex_to_tree.cxx
 *  This file implements the routines for the t8_cmesh_conn_vertex_to_tree_c struct.
 */

t8_cmesh_tree_vertex_list &
t8_cmesh_vertex_conn_vertex_to_tree_c::get_tree_list_of_vertex (t8_gloidx_t global_vertex_id)
{
  T8_ASSERT (0 <= global_vertex_id);

  /* Use at() to look for the vertex entry.
   * If the entry does not exist an exception of
   * type std::out_of_range is thrown. */

  try {
    return vertex_to_tree.at (global_vertex_id);
  } catch (const std::out_of_range &e) {
    t8_errorf ("ERROR: Could not find vertex %li for cmesh.\n", global_vertex_id);
    SC_ABORTF ("Caught exception 'out of range': %s\n", e.what ());
  }
}

void
t8_cmesh_vertex_conn_vertex_to_tree_c::add_vertex_to_tree (t8_cmesh_t cmesh, t8_gloidx_t global_vertex_id,
                                                           t8_locidx_t ltreeid, int tree_vertex)
{
  T8_ASSERT (0 <= global_vertex_id);
  T8_ASSERT (t8_cmesh_treeid_is_local_tree (cmesh, ltreeid) || t8_cmesh_treeid_is_ghost (cmesh, ltreeid));
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

#if T8_ENABLE_DEBUG
  t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, ltreeid);
  int num_tree_vertices = t8_eclass_num_vertices[tree_class];
  T8_ASSERT (0 <= tree_vertex && tree_vertex < num_tree_vertices);
#endif

  t8_cmesh_tree_vertex_pair pair (ltreeid, tree_vertex);
  t8_cmesh_tree_vertex_list &list_of_globalid = vertex_to_tree[global_vertex_id];

  list_of_globalid.push_back (pair);
}