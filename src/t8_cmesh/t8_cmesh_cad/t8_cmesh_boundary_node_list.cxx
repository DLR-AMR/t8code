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

/** \file t8_cmesh_boundary_node_list.cxx
 *
 * TODO: document this file
 */

#include <t8_cmesh.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <span>

t8_boundary_node_list::t8_boundary_node_list (t8_cmesh_t cmesh_in): cmesh (cmesh_in)
{
  boundary_node_list = compute_boundary_node ();
}

std::unordered_set<t8_gloidx_t>
t8_boundary_node_list::compute_boundary_node ()
{
  std::unordered_set<t8_gloidx_t> boundary_node_list = {};

  const t8_locidx_t num_trees = cmesh->num_local_trees;

  /* Iterate through trees */
  for (t8_locidx_t i_tree = 0; i_tree < num_trees; i_tree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, i_tree);
    int num_faces = t8_eclass_num_faces[(int) eclass];
    const std::span<const t8_gloidx_t> global_vertices_of_tree
      = cmesh->vertex_connectivity->get_global_vertices_of_tree (i_tree); /* Get global node IDs of i_tree */
    for (int i_face = 0; i_face < num_faces; i_face++) {                  /* Iterate through faces of i_tree */
      const int vertex_per_face = t8_eclass_num_vertices[t8_eclass_face_types[eclass][i_face]];

      if (t8_cmesh_tree_face_is_boundary (cmesh, i_tree, i_face)) { /* Check if i_face is boundary face*/
        for (int count = 0; count < vertex_per_face; count++) {
          boundary_node_list.insert (
            global_vertices_of_tree[t8_face_vertex_to_tree_vertex[eclass][i_face]
                                                                 [count]]); /* Append global node IDs of i_face*/
        }
      }
    }
  }
  return boundary_node_list;
}

std::unordered_set<t8_gloidx_t>
t8_boundary_node_list::get_boundary_node_list ()
{
  return this->boundary_node_list;
}
