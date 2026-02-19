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

/** \file t8_cmesh_edge_connectivity.hxx
 * We define classes and interfaces for a global edge enumeration
 * of a cmesh.
 */

#ifndef T8_CMESH_EDGE_CONNECTIVITY
#define T8_CMESH_EDGE_CONNECTIVITY

#include <memory>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_conn_edge_to_tree.hxx>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_conn_tree_to_edge.hxx>

/**
 * A class to hold the edge connectivity of a cmesh.
 */
struct t8_cmesh_edge_connectivity
{
 public:
  /**
   * Constructor.
   */
  t8_cmesh_edge_connectivity ()
    : current_state (state::INITIALIZED), global_number_of_edges (0), local_number_of_edges (0),
      associated_cmesh (nullptr) {};

  /**
   * Destructor.
   */
  ~t8_cmesh_edge_connectivity () {};

  /** The state this connectivity can be in. */
  enum struct state {
    INITIALIZED,          /**< Initialized but not filled */
    TREE_TO_EDGE_VALID, /**< Ready to use, but only tree_to_edge functionality. */
    VALID                 /**< Ready to use for full edge connectivity. Cannot be altered anymore. */
  };

  /* Setter functions */

  /** Set all global edge ids of a local tree.
    * \param [in] cmesh The considered cmesh
    * \param [in] global_tree A global tree id of \a cmesh
    * \param [in] global_tree_edges The ids of the global edges in order of \a local_tree's edges.
    * \param [in] num_edges Must match the number of edges of \a local_tree
    *
    * \note \a cmesh must not be committed.
    */
  inline void
  set_global_edge_ids_of_tree_edges (const t8_cmesh_t cmesh, const t8_gloidx_t global_tree,
                                          const t8_gloidx_t *global_tree_edges, const int num_edges)
  {
    T8_ASSERT (t8_cmesh_is_initialized (cmesh));
    T8_ASSERT (!t8_cmesh_is_committed (cmesh));
    T8_ASSERT (current_state != state::VALID);
    current_state = state::TREE_TO_EDGE_VALID;
    if (associated_cmesh == nullptr) {
      associated_cmesh = cmesh;
    }
    T8_ASSERT (associated_cmesh == cmesh);
    return tree_to_edge.set_global_edge_ids_of_tree_edges (cmesh, global_tree, global_tree_edges,
                                                                  num_edges);
  }

  /** Function to fill ett from the stored tte information.
   * Sets all global ids and associated tree edges from
   * the associated cmesh.
   * Afterwards, this class is ready to be used and cannot be altered.
  */
  void
  build_edge_to_tree ()
  {
    T8_ASSERT (current_state == state::TREE_TO_EDGE_VALID);
    edge_to_tree.build_from_tte (associated_cmesh, tree_to_edge);
    global_number_of_edges = edge_to_tree.edge_to_tree.size ();
    current_state = state::VALID;
  }

  /* Getter functions */

  /** Return the total number of global edges of a cmesh (across all processes).
   * \return The total number of global edges of \a cmesh.
  */
  inline t8_gloidx_t
  get_global_number_of_edges () const
  {
    return global_number_of_edges;
  }

  /** Return the number of process local global edges of a cmesh.
   * \return The number of process local global edges of \a cmesh.
  */
  inline t8_gloidx_t
  get_local_number_of_edges ()
  {
    /* TODO: Change this after the global edge connectivity is working partitioned. */
    return local_number_of_edges = global_number_of_edges;
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

  /** Given a global edge id, return the list of trees that share this edge.
   * \param [in] edge_id A global edge id.
   * \return The trees and their local edge ids matching \a edge_id.
   */
  inline const tree_edge_list &
  edge_to_trees (const t8_gloidx_t edge_id) const
  {
    T8_ASSERT (edge_to_tree.is_committed ());
    T8_ASSERT (current_state == state::VALID);
    return edge_to_tree.get_tree_list_of_edge (edge_id);
  }

  /** Given a local tree (or ghost) return the list of global edges
   * this tree has (in order of its local edges).
   * \param [in] ltree A local tree id.
   * \return The list of global edges that this tree has.
   */
  inline const std::span<const t8_gloidx_t>
  tree_to_edges (const t8_locidx_t ltree) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_EDGE_VALID);
    T8_ASSERT (tree_to_edge.get_state () == t8_cmesh_edge_conn_tree_to_edge::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_edge.get_global_edges (associated_cmesh, ltree);
  }

  /** Given a local tree (or ghost) and a local edge of that tree
   * return the global edge id of that edge.
   * \param [in] ltree A local tree.
   * \param [in] ltree_edge A local edge of \a ltree.
   * \return The global edge associated with \a ltree_edge.
   */
  inline t8_gloidx_t
  treeedge_to_edge (const t8_locidx_t ltree, const t8_locidx_t ltree_edge) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_EDGE_VALID);
    T8_ASSERT (tree_to_edge.get_state () == t8_cmesh_edge_conn_tree_to_edge::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_edge.get_global_edge (associated_cmesh, ltree, ltree_edge);
  }

  /** Get the global edge indices of a tree in its local edge order.
   * \param [in] local_tree A local tree in \a cmesh.
   * \return The global edges of \a local_tree
   */
  inline const std::span<const t8_gloidx_t>
  get_global_edges_of_tree (const t8_locidx_t local_tree) const
  {
    T8_ASSERT (current_state >= state::TREE_TO_EDGE_VALID);
    T8_ASSERT (tree_to_edge.get_state () == t8_cmesh_edge_conn_tree_to_edge::state::FILLED);
    T8_ASSERT (associated_cmesh != nullptr);
    return tree_to_edge.get_global_edges (associated_cmesh, local_tree);
  }

  /** Get a single global edge index of a local tree's local edge.
  * \param [in] local_tree A local tree in \a cmesh.
  * \param [in] local_tree_edge A local edge of \a local_tree
  * \return The global edge matching \a local_tree_edge of \a local_tree.
  */
  inline t8_gloidx_t
  get_global_edge_of_tree (const t8_locidx_t local_tree, const int local_tree_edge) const
  {
    T8_ASSERT (associated_cmesh != nullptr);
    auto edges_of_tree = get_global_edges_of_tree (local_tree);
    T8_ASSERT (edges_of_tree.size () > static_cast<size_t> (local_tree_edge));
    return edges_of_tree[local_tree_edge];
  }

  /** Get the list of global trees and local edge ids a global edge is connected to.
   *
   * \param [in] global_edge_id The global id of a edge in the cmesh.
   * \return The list of global tree ids and local edge ids of \a global_edge_id.
   */
  inline const tree_edge_list &
  get_tree_list_of_edge (const t8_gloidx_t global_edge_id) const
  {
    T8_ASSERT (current_state == state::VALID);
    return edge_to_tree.get_tree_list_of_edge (global_edge_id);
  }

  /** Get the number of global trees a global edge is connected to.
   * \param [in] global_edge_id The global id of a edge in the cmesh.
   * \note if a tree is contained multiple times it is counted as multiple entries.
   * Example: For a quad where all 4 edges map to a single global edge this function will return 4.
   */
  inline int
  get_num_trees_at_edge (const t8_gloidx_t global_edge_id) const
  {
    T8_ASSERT (current_state == state::VALID);
    return get_tree_list_of_edge (global_edge_id).size ();
  }

 private:
  /** The internal state. Indicating whether this structure is new and unfilled, the tte was filled or the edge conn was built completely. */
  state current_state;

  /** The process global number of global edges. */
  t8_gloidx_t global_number_of_edges;

  /** Currently not used/equal to global number of edges */
  /* TODO: This needs to be properly implemented when the global edge list works for partitioned cmeshes. */
  t8_locidx_t local_number_of_edges;

  /** The internal edge to tree storage. */
  t8_cmesh_edge_conn_edge_to_tree edge_to_tree;

  /** The internal tree to edge storage. */
  t8_cmesh_edge_conn_tree_to_edge tree_to_edge;

  /** A pointer to the cmesh for attribute retrieval */
  t8_cmesh_t associated_cmesh;
};

#endif /* !T8_CMESH_EDGE_CONNECTIVITY */
