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

/** \file t8_cmesh_vertex_conn_tree_to_vertex.hxx
 * Class to save data structure for tree_to_vertex_lists of the cmesh.
 * When the cmesh stores global vertex numbers, we require a lookup that
 * matches a tree and its local vertex to a global vertex id.
 * This lookup is encoded in the t8_cmesh_vertex_conn_tree_to_vertex struct.
 */

/* TODO:
 *  It is probably best to set all global ids of a single tree as one attribute.
 * That way we can store it as a single arrays of id's.
 * We will probably most often want to access all ids of one tree, so a function returning an array of ids
 * is a good idea anyway and if we already store them as such then we do not need to do any data movement
 * when accessing.
 *
 * On the downside we will only have a "set all ids of a tree" function and no "set this single id for this tree and this vertex" function.
 */

#ifndef T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX
#define T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX

#include <algorithm>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_conn_vertex_to_tree.hxx>

/* forward declaration of ttv class needed since the two class headers include each other. */
class t8_cmesh_vertex_conn_vertex_to_tree;

class t8_cmesh_vertex_conn_tree_to_vertex {
 public:
  /** Standard constructor. Does nothing. */
  t8_cmesh_vertex_conn_tree_to_vertex (): state (EMPTY)
  {
  }

  /** Constructor from a cmesh where all the attributes are set.
   * Currently unclear if we implement this eventually.
   * If we do so: Should the cmesh be already committed, or in pre-commit state but attributes set?
   * 
   * \note This function is not implemented yet.
   */
  t8_cmesh_vertex_conn_tree_to_vertex (const t8_cmesh_t cmesh)
  {
    SC_ABORT ("not implemented.");
  }

  /** Constructor from a cmesh and a given vertex to tree connectivity.
   *
   * \param [in] cmesh_from A committed cmesh.
   * \param [in] cmesh      An initialized but not committed cmesh that is to be derived from \a cmesh_from.
   * \param [in] vtt        A committed vertex to tree connectivity for \a cmesh_from.
   * 
   * As a result a tree to vertec connectivity for \a cmesh will be constructed.
   * \note \a cmesh_from must be committed.
   * \note \a cmesh must not be committed.
   * \note \a vtt must be committed.
   * \note This does not work until issue #923 https://github.com/DLR-AMR/t8code/issues/923 is resolved.
   */
  t8_cmesh_vertex_conn_tree_to_vertex (const t8_cmesh_t cmesh_from, const t8_cmesh_t cmesh,
                                       const struct t8_cmesh_vertex_conn_vertex_to_tree &vtt);

  /* Setter functions */
  /** Set all global vertex ids of a local tree.
   * \param [in] cmesh The considered cmesh
   * \param [in] local_tree A local tree id of \a cmesh
   * \param [in] global_vertex_id The ids of the global vertices in order of \a local_tree's vertices.
   * \param [in] num_vertices Must match the number of vertices of \a local_tree
   *
   * \note \a cmesh must not be committed.
  */
  void
  set_global_vertex_ids_of_tree_vertices (const t8_cmesh_t, const t8_gloidx_t global_tree,
                                          const t8_gloidx_t *global_tree_vertices, const int num_vertices);

  t8_gloidx_t
  get_global_vertex (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_vertex,
                     const int num_tree_vertices) const;

  const t8_gloidx_t *
  get_global_vertices (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int num_vertices) const;

  const int
  get_state ()
  {
    return state;
  }

  friend struct t8_cmesh_vertex_connectivity;

 private:
  enum {
    EMPTY, /*< Is initialized but empty. */
    FILLED /*< Is filled with at least one entry. */
  } state;
};

#endif /* !T8_CMESH_VERTEX_CONN_TREE_TO_VERTEX_HXX */
