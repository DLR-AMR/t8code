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

/** \file t8_cmesh_vertex_connectivity.hxx
 * We define classes and interfaces for a global vertex enumeration
 * of a cmesh.
 */

#ifndef T8_CMESH_VERTEX_CONNECTIVITY_HXX
#define T8_CMESH_VERTEX_CONNECTIVITY_HXX

#include <t8_cmesh.h>

typedef struct t8_cmesh_vertex_connectivity
{
 public:
  /* Setter functions */
  void
  set_state (const t8_cmesh_vertex_connectivity_state_t state);

  /* Given a cmesh, build up the vertex_to_tree and tree_to_vertex members.
   * \return: some error value to be specified.
   * On error, \state will be set to ERROR. 
   * The cmesh must not be committed, but all tree information and neighbor information must
   * have been set. 
   * Currently, \a cmesh has to be replicated. */
  void
  build (const t8_cmesh_t cmesh);

  const t8_gloidx_t
  get_global_number_of_vertices ();

  const t8_gloidx_t
  get_local_number_of_vertices ();

  /* Getter functions */

  /* Return the state of the connectivity */
  const t8_cmesh_vertex_connectivity_state_t
  get_state ();

  /* Given a global vertex id, return the list of trees that share this vertex */
  const *DEFINE_DATA_TYPE
  vertex_to_trees (t8_gloidx_t vertex_id);

  /* Given a local tree (or ghost) return the list of global vertices
   * this tree has (in order of its local vertices). */
  const *ARRAY_OF_GLOIDX_T
  tree_to_vertices (t8_locidx_t ltree);

  /* Given a local tree (or ghost) and a local vertex of that tree
   * return the global vertex id of that vertex */
  const t8_gloidx_t
  treevertex_to_vertex (t8_locidx_t ltree, t8_locidx_t ltree_vertex);

 private:
  t8_cmesh_vertex_connectivity_state_t state;

  t8_gloidx_t global_number_of_vertices;

  /* Currently not used/equal to global number of vertices */
  t8_locidx_t local_number_of_vertices;

  t8_cmesh_vertex_conn_vertex_to_tree_c vertex_to_tree;

  t8_cmesh_vertex_conn_tree_to_vertex_c tree_to_vertex;

} t8_cmesh_vertex_connectivity_c;

#endif /* !T8_CMESH_VERTEX_CONNECTIVITY_HXX */
