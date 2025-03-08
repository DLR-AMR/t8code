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

/** \file t8_cmesh_vertex_connectivity.cxx
 * We define classes and interfaces for a global vertex enumeration
 * of a cmesh.
 */

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity.hxx>

/*
Use case

initialize cmeesh
Add trees
Add tree to vertex information
Start commit cmesh
  commit cmesh internally
  build vertex to tree information
Finish commit cmesh

Access vtt and ttv information

-> ttv must be added before commit
global tree id -> vertex_list

-> vtt must be added after commit (not by user)

*/

/* Setter functions */

void
t8_cmesh_set_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                      const t8_gloidx_t *global_tree_vertices, const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  cmesh->vertex_connectivity->set_global_vertex_ids_of_tree_vertices (cmesh, global_tree, global_tree_vertices,
                                                                      num_vertices);
}

t8_gloidx_t
t8_cmesh_get_num_global_vertices (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  return cmesh->vertex_connectivity->get_global_number_of_vertices ();
}

t8_locidx_t
t8_cmesh_get_num_local_vertices (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  return cmesh->vertex_connectivity->get_local_number_of_vertices ();
}

const t8_gloidx_t *
t8_cmesh_get_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  return cmesh->vertex_connectivity->get_global_vertices_of_tree (cmesh, local_tree, num_vertices);
}

const t8_gloidx_t
t8_cmesh_get_global_vertex_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_vertex,
                                    const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const t8_gloidx_t *vertices_of_tree = t8_cmesh_get_global_vertices_of_tree (cmesh, local_tree, num_vertices);
  return vertices_of_tree[local_tree_vertex];
}

const t8_cmesh_vertex_conn_vertex_to_tree::tree_vertex_list &
t8_cmesh_get_vertex_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_vertex)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->vertex_connectivity->get_tree_list_of_vertex (global_vertex);
}

const int
t8_cmesh_get_num_trees_at_vertex (const t8_cmesh_t cmesh, t8_gloidx_t global_vertex)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->vertex_connectivity->get_tree_list_of_vertex (global_vertex).size ();
}
