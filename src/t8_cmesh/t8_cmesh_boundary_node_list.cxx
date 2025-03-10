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
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh_boundary_node_list.hxx>
#include <unordered_set>

t8_boundary_node_list::t8_boundary_node_list(t8_cmesh_t cmesh) : cmesh(cmesh){
  boundary_node_list = compute_boundary_node(cmesh);
}

std::unordered_set<t8_gloidx_t> t8_boundary_node_list::compute_boundary_node(t8_cmesh_t cmesh){
    std::unordered_set<t8_gloidx_t> boundary_node_list = {};
    const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees(cmesh);
    
    for (t8_locidx_t i_tree = 0; i_tree < num_trees; i_tree++){
      int8_t *ttf;
      t8_locidx_t *face_neighbor;
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, i_tree);
      int num_faces = t8_eclass_num_faces[(int) eclass];
      (void) t8_cmesh_trees_get_tree_ext (cmesh->trees, i_tree, &face_neighbor, &ttf);

      for (int i_face = 0; i_face < num_faces; i_face++)
        if (face_neighbor[i_face] == i_tree && ttf[i_face] == i_face) {
        /* The tree is connected to itself at the same face.
        * Thus this is a domain boundary */
          boundary_node_list.insert(t8_cmesh_get_global_id(cmesh, i_tree));
      }
  }
    return boundary_node_list;
}

std::unordered_set<t8_gloidx_t> t8_boundary_node_list::get_boundary_node_list() {
    return this->boundary_node_list;
}