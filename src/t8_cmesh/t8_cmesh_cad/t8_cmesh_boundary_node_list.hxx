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

/** \file t8_cmesh_boundary_node_list.hxx
 * Defines the boundary node list, which computes the boundary nodes of a cmesh.
 */

#ifndef T8_CMESH_BOUNDARY_NODE_LIST_HXX
#define T8_CMESH_BOUNDARY_NODE_LIST_HXX

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <unordered_set>

class t8_boundary_node_list {
 public:
  t8_boundary_node_list (t8_cmesh_t cmesh_in);

  std::unordered_set<t8_gloidx_t>
  get_boundary_node_list ();

 private:
  std::unordered_set<t8_gloidx_t>
  compute_boundary_node ();

  t8_cmesh_t cmesh;
  std::unordered_set<t8_gloidx_t> boundary_node_list;
};

#endif /* T8_CMESH_BOUNDARY_NODE_LIST_HXX */
