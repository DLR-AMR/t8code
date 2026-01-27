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

/** \file t8_cmesh_vertex_connectivity.hxx
 * We define classes and interfaces for a global vertex enumeration
 * of a cmesh.
 */

#ifndef T8_CMESH_VERTEX_CONNECTIVITY
#define T8_CMESH_VERTEX_CONNECTIVITY

#include <memory>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.hxx>

/**
 * A class to hold the vertex connectivity of a cmesh.
 */
struct t8_cmesh_vertex_connectivity
{
 public:
  /**
   * Constructor.
   */
  t8_cmesh_vertex_connectivity ()
    : current_state (state::INITIALIZED), global_number_of_vertices (0), local_number_of_vertices (0),
      associated_cmesh (nullptr) {};

  /**
   * Destructor.
   */
  ~t8_cmesh_vertex_connectivity () {};

  /** The state this connectivity can be in. */
  enum struct state {
    INITIALIZED,          /**< Initialized but not filled */
    TREE_TO_VERTEX_VALID, /**< Ready to use, but only tree_to_vertex functionality. */
    VALID                 /**< Ready to use for full vertex connectivity. Cannot be altered anymore. */
  };

  /* Setter functions */

  /** Set all global vertex ids of a local tree.
    * \param [in] cmesh The considered cmesh
    * \param [in] global_tree A global tree id of \a cmesh
    * \param [in] global_tree_vertices The ids of the global vertices in order of \a local_tree's vertices.
    * \param [in] num_vertices Must match the number of vertices of \a local_tree
    *
    * \note \a cmesh must not be committed.
    */
  inline void
  set_global_vertex_ids_of_tree_vertices (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                          const t8_gloidx_t *global_tree_vertices, const int num_vertices)
  {
    T8_ASSERT (t8_cmesh_is_initialized (cmesh));
    T8_ASSERT (!t8_cmesh_is_committed (cmesh));
    T8_ASSERT (current_state != state::VALID);
    current_state = state::TREE_TO_VERTEX_VALID;
    if (associated_cmesh == nullptr) {
      associated_cmesh = cmesh;
    }
    T8_ASSERT (associated_cmesh == cmesh);
    return tree_to_vertex.set_global_vertex_ids_of_tree_vertices (cmesh, global_tree, global_tree_vertices,
                                                                  num_vertices);
  }

  /** Function to fill vtt from the stored ttv information.
   * Sets all global ids and associated tree vertices from
   * the associated cmesh.
   * Afterwards, this class is ready to be used and cannot be altered.
  */
  void
  build_vertex_to_tree ()
  {
    T8_ASSERT (current_state == state::TREE_TO_VERTEX_VALID);
    vertex_to_tree.build_from_ttv (associated_cmesh, tree_to_vertex);
    global_number_of_vertices = vertex_to_tree.vertex_to_tree.size ();
    current_state = state::VALID;
  }

  /* Getter functions */

  /** Return the total number of global vertices of a cmesh (across all processes).
   * \return The total number of global vertices of \a cmesh.
  */
  inline t8_gloidx_t
  get_global_number_of_vertices () const
  {
    return global_number_of_vertices;
  }

  /** Return the number of process local global vertices of a cmesh.
   * \return The number of process local global vertices of \a cmesh.
  */
  inline t8_gloidx_t
  get_local_number_of_vertices ()
  {
    /* TODO: Change this after the global vertex connectivity is working partitioned. */
    return local_number_of_vertices = global_number_of_vertices;
  }

  /** Return the state of the connectivity.
   * \return The current state (initialized, committed, etc) of this object.
   * \ref state
   */
  inline state
  get_state () const
  {
    return current_state;
  }

  /** Given a global vertex id, return the list of trees that share this vertex.
   * \param [in] vertex_id A global vertex id.
   * \return The trees and their local vertex ids matching \a vertex_id.
   */
  inline const tree_vertex_list &
  vertex_to_trees (const t8_gloidx_t vertex_id) const
  {
    T8_ASSERT (vertex_to_tree.is_committed ());
    T8_ASSERT (current_state == state::VALID);
    return vertex_to_tree.get_tree_list_of_vertex (vertex_id);
  }

  /** Given a local tree (or ghost) return the list of global vertices
   * this tree has (in order of its local vertices).
   * \param [in] ltree A local tree id.
   * \return The list of global vertices that this tree has.
   */
  inline const std::span<const t8_gloidx_t>
  tree_to_vertices (const t8_locidx_t ltree) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_VERTEX_VALID);
    T8_ASSERT (tree_to_vertex.get_state () == t8_cmesh_vertex_conn_tree_to_vertex::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_vertex.get_global_vertices (associated_cmesh, ltree);
  }

  /** Given a local tree (or ghost) and a local vertex of that tree
   * return the global vertex id of that vertex.
   * \param [in] ltree A local tree.
   * \param [in] ltree_vertex A local vertex of \a ltree.
   * \return The global vertex associated with \a ltree_vertex.
   */
  inline t8_gloidx_t
  treevertex_to_vertex (const t8_locidx_t ltree, const t8_locidx_t ltree_vertex) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_VERTEX_VALID);
    T8_ASSERT (tree_to_vertex.get_state () == t8_cmesh_vertex_conn_tree_to_vertex::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_vertex.get_global_vertex (associated_cmesh, ltree, ltree_vertex);
  }

  /** Get the global vertex indices of a tree in its local vertex order.
   * \param [in] local_tree A local tree in \a cmesh.
   * \return The global vertices of \a local_tree
   */
  inline const std::span<const t8_gloidx_t>
  get_global_vertices_of_tree (const t8_locidx_t local_tree) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_VERTEX_VALID);
    T8_ASSERT (tree_to_vertex.get_state () == t8_cmesh_vertex_conn_tree_to_vertex::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_vertex.get_global_vertices (associated_cmesh, local_tree);
  }

  /** Get a single global vertex index of a local tree's local vertex.
  * \param [in] local_tree A local tree in \a cmesh.
  * \param [in] local_tree_vertex A local vertex of \a local_tree
  * \return The global vertex matching \a local_tree_vertex of \a local_tree.
  */
  inline t8_gloidx_t
  get_global_vertex_of_tree (const t8_locidx_t local_tree, const int local_tree_vertex) const
  {
    T8_ASSERT (associated_cmesh != nullptr);
    auto vertices_of_tree = get_global_vertices_of_tree (local_tree);
    T8_ASSERT (vertices_of_tree.size () > static_cast<size_t> (local_tree_vertex));
    return vertices_of_tree[local_tree_vertex];
  }

  /** Get the list of global trees and local vertex ids a global vertex is connected to.
   *
   * \param [in] global_vertex_id The global id of a vertex in the cmesh.
   * \return The list of global tree ids and local vertex ids of \a global_vertex_id.
   */
  inline const tree_vertex_list &
  get_tree_list_of_vertex (const t8_gloidx_t global_vertex_id) const
  {
    T8_ASSERT (current_state == state::VALID);
    return vertex_to_tree.get_tree_list_of_vertex (global_vertex_id);
  }

  /** Get the number of global trees a global vertex is connected to.
   * \param [in] global_vertex_id The global id of a vertex in the cmesh.
   * \note if a tree is contained multiple times it is counted as multiple entries.
   * Example: For a quad where all 4 vertices map to a single global vertex this function will return 4.
   */
  inline int
  get_num_trees_at_vertex (const t8_gloidx_t global_vertex_id) const
  {
    T8_ASSERT (current_state == state::VALID);
    return get_tree_list_of_vertex (global_vertex_id).size ();
  }

 private:
  /** The internal state. Indicating whether this structure is new and unfilled, the ttv was filled or the vertex conn was built completely. */
  state current_state;

  /** The process global number of global vertices. */
  t8_gloidx_t global_number_of_vertices;

  /** Currently not used/equal to global number of vertices */
  /* TODO: This needs to be properly implemented when the global vertex list works for partitioned cmeshes. */
  t8_locidx_t local_number_of_vertices;

  /** The internal vertex to tree storage. */
  t8_cmesh_vertex_conn_vertex_to_tree vertex_to_tree;

  /** The internal tree to vertex storage. */
  t8_cmesh_vertex_conn_tree_to_vertex tree_to_vertex;

  /** A pointer to the cmesh for attribute retrieval */
  t8_cmesh_t associated_cmesh;
};

#endif /* !T8_CMESH_VERTEX_CONNECTIVITY */
