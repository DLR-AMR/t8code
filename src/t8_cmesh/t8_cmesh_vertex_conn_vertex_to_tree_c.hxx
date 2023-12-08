/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_vertex_conn_vertex_to_tree_c.hxx
 * Class to save data structure for vertex_to_tree_lists
 */

#ifndef T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX
#define T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX

#include <t8_cmesh.h>

typedef struct t8_cmesh_vertex_conn_vertex_to_tree_c
{
 public:
  /* Setter functions */
  /* Given a cmesh, build up the vertex_to_tree.
   * \return: some error value to be specified.
   * The cmesh must not be committed, but all tree information and neighbor information must
   * have been set. 
   * Currently, \a cmesh has to be replicated. */
  void
  set_vertex_to_tree_list ( const t8_cmesh_t cmesh );

  vector< t8_locidx_t >
  get_tree_list_of_vertex( t8_gloidx_t vertex_id );

 private:
  /* Vector of vectors: For each vertex one list of trees */
  vector<vector< t8_locidx_t >> vertex_to_tree_list;

  /* Setter functions */
  /* A single value is added to the vertex_to_tree_list */
  void
  set_value_vertex_to_tree_list ( t8_gloidx_t vertex_id, t8_locidx_t treeid );


} T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX;

#endif /* !T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX */
