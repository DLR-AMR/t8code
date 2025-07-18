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
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>

/* Setter functions */

void
t8_cmesh_set_vertex_conn (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  if (cmesh->vertex_connectivity == nullptr)
    cmesh->vertex_connectivity = new t8_cmesh_vertex_connectivity ();
}

int
t8_cmesh_get_vertex_conn_status (t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  if (cmesh->vertex_connectivity == nullptr)
    return 0;
  if (cmesh->committed) {
    return cmesh->vertex_connectivity->get_state () == t8_cmesh_vertex_connectivity::state::VALID;
  }
  return 1;
}

void
t8_cmesh_set_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                      const t8_gloidx_t *global_tree_vertices, const int num_vertices)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!t8_cmesh_is_committed (cmesh));
  if (cmesh->vertex_connectivity == nullptr)
    cmesh->vertex_connectivity = new t8_cmesh_vertex_connectivity ();
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
  T8_ASSERTF (!t8_cmesh_is_partitioned (cmesh), "Global vertex ids currently not supported for partitioned cmeshes.");
  return cmesh->vertex_connectivity->get_local_number_of_vertices ();
}

const t8_gloidx_t *
t8_cmesh_get_global_vertices_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, int *num_vertices)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const auto global_vertices = cmesh->vertex_connectivity->get_global_vertices_of_tree (local_tree);
  *num_vertices = global_vertices.size ();
  return global_vertices.data ();
}

t8_gloidx_t
t8_cmesh_get_global_vertex_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_vertex)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  int num_vertices;
  const t8_gloidx_t *vertices_of_tree = t8_cmesh_get_global_vertices_of_tree (cmesh, local_tree, &num_vertices);
  T8_ASSERT (num_vertices > local_tree_vertex);
  return vertices_of_tree[local_tree_vertex];
}

const tree_vertex_list &
t8_cmesh_get_vertex_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_vertex)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->vertex_connectivity->get_tree_list_of_vertex (global_vertex);
}

int
t8_cmesh_get_num_trees_at_vertex (const t8_cmesh_t cmesh, t8_gloidx_t global_vertex)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->vertex_connectivity->get_tree_list_of_vertex (global_vertex).size ();
}
