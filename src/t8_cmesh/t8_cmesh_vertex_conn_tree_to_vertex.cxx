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

#include <t8_cmesh_conn_tree_to_vertex.hxx>

/** \file t8_cmesh_conn_tree_to_vertex.cxx
 *  This file implements the routines for the t8_cmesh_conn_tree_to_vertex_c struct.
 */

  /* Set all global vertex ids of a local tree. 
   * \param[in] cmesh The considered cmesh
   * \param[in] local_tree A local tree id of \a cmesh
   * \param[in] global_vertex_id The ids of the global vertices in order of \a local_tree's vertices.
   * \param[in] num_vertices Must match the number of vertices of \a local_tree
   * 
   * \note Cmesh must not be committed.
  */
  void
  set_global_vertex_ids_of_tree_vertices (const t8_cmesh_t, const t8_locidx_t local_tree, const int *global_tree_vertex, const int num_vertices)
  {

  }