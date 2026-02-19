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

/** \file t8_cmesh_edge_connectivity.cxx
 * We define classes and interfaces for a global edge enumeration
 * of a cmesh.
 */

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_connectivity.hxx>

/* Setter functions */

void
t8_cmesh_set_global_edges_of_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                      const t8_gloidx_t *global_tree_edges, const int num_edges)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  cmesh->edge_connectivity->set_global_edge_ids_of_tree_edges (cmesh, global_tree, global_tree_edges,
                                                                      num_edges);
}

t8_gloidx_t
t8_cmesh_get_num_global_edges (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  return cmesh->edge_connectivity->get_global_number_of_edges ();
}

t8_locidx_t
t8_cmesh_get_num_local_edges (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERTF (!t8_cmesh_is_partitioned (cmesh), "Global edge ids currently not supported for partitioned cmeshes.");
  return cmesh->edge_connectivity->get_local_number_of_edges ();
}

const t8_gloidx_t *
t8_cmesh_get_global_edges_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, int *num_edges)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  const auto global_edges = cmesh->edge_connectivity->get_global_edges_of_tree (local_tree);
  if (num_edges) {
    *num_edges = global_edges.size ();
  }
  return global_edges.data ();
}

t8_gloidx_t
t8_cmesh_get_global_edge_of_tree (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_edge)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  int num_edges;
  const t8_gloidx_t *edges_of_tree = t8_cmesh_get_global_edges_of_tree (cmesh, local_tree, &num_edges);
  T8_ASSERT (num_edges > local_tree_edge);
  return edges_of_tree[local_tree_edge];
}

const tree_edge_list &
t8_cmesh_get_edge_to_tree_list (const t8_cmesh_t cmesh, const t8_gloidx_t global_edge)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->edge_connectivity->get_tree_list_of_edge (global_edge);
}

int
t8_cmesh_get_num_trees_at_edge (const t8_cmesh_t cmesh, t8_gloidx_t global_edge)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));

  return cmesh->edge_connectivity->get_tree_list_of_edge (global_edge).size ();
}

int
t8_cmesh_uses_edge_connectivity (const t8_cmesh_t cmesh)
{
  return (cmesh->edge_connectivity->get_state () == t8_cmesh_edge_connectivity::state::TREE_TO_EDGE_VALID)
         || (cmesh->edge_connectivity->get_state () == t8_cmesh_edge_connectivity::state::VALID);
}
