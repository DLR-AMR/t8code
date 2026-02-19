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

/** \file t8_cmesh_edge_conn_edge_to_tree.hxx
 * Class to save data structures cmesh for edge_to_tree_lists
 */

#ifndef T8_CMESH_EDGE_CONN_EDGE_TO_TREE_HXX
#define T8_CMESH_EDGE_CONN_EDGE_TO_TREE_HXX

#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_conn_tree_to_edge.hxx>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_connectivity_types.hxx>

/** forward declaration of tte class needed since the two class headers include each other. */
struct t8_cmesh_edge_conn_tree_to_edge;

/** This class stores the edge to tree lookup for
 * global edge indices for a cmesh.
 * Thus, given a global edge id the class provides
 * information about the trees the edge belongs to and
 * the corresponding local edge ids inside these trees.
 *
 * In particular, this class stores the lookup
 *
 * global_edge_id -> List of (tree, tree_local_edge)
 *
 * for a cmesh.
 * It is the opposite lookup as \ref t8_cmesh_edge_conn_tree_to_edge.hxx
 *
 * The global edge ids must not be contiguous, that is, we have some set
 *
 * {I_0 < I_1 < ...< I_N} of natural numbers corresponding to the N+1 edges.
 *
 * I_0 does not have to be 0 and I_N does not have to be N.
 *
 * So we need lookup: I_i -> i
 *  store this in a hash table.
 *
 * datatypes:
 *
 * global id: t8_gloidx_t
 * (tree_id, tree_edge): std::pair<t8_locidx_t, int> = TV_PAIR
 * List of (tree_id, tree_edge): std::vector<PAIR> = TV_LIST
 * Table global_id -> TV_LIST: std::unordered_map<t8_gloidx_t, TV_LIST>
 *
 */
struct t8_cmesh_edge_conn_edge_to_tree
{
 public:
  /** Standard constructor.
   * Initializes the class and allows setting edge entries
   * via \ref add_edge_to_tree
   */
  t8_cmesh_edge_conn_edge_to_tree (): state (INITIALIZED)
  {
  }

  /** Function to fill ett from cmesh and tte information.
   * Sets all global ids and associated tree edges from
   * the given input class.
   * Afterwards, the class is set to committed and can be used.
   *
   * \param [in] cmesh A committed cmesh with set tree to edge entries.
   * \param [in] tte A filled tree to edge list for \a cmesh.
  */
  void
  build_from_tte (const t8_cmesh_t cmesh, t8_cmesh_edge_conn_tree_to_edge& tte);

  /* Setter functions */
  /** Given a cmesh, build up the edge_to_tree.
   * \param [in] cmesh An initialized but not yet committed cmesh.
   * The cmesh must not be committed, but all tree information and neighbor information must
   * have been set.
   * Currently, \a cmesh has to be replicated. */
  void
  set_edge_to_tree_list (const t8_cmesh_t cmesh);

  /** Get the list of global trees and local edge ids a global edge is connected to.
   *
   * \param [in] global_edge_id The global id of a edge in the cmesh.
   * \return The list of local tree ids and local edge ids of \a global_edge_id.
   */
  inline const tree_edge_list&
  get_tree_list_of_edge (const t8_gloidx_t global_edge_id) const
  {
    T8_ASSERT (is_committed ());
    T8_ASSERT (0 <= global_edge_id);

    /* Use at() to look for the edge entry.
    * If the entry does not exist an exception of
    * type std::out_of_range is thrown. */

    try {
      return edge_to_tree.at (global_edge_id);
    } catch (const std::out_of_range& e) {
      t8_errorf ("ERROR: Could not find edge %" T8_GLOIDX_FORMAT " for cmesh.\n", global_edge_id);
      SC_ABORTF ("Caught exception 'out of range': %s\n", e.what ());
    }
  }

  /** Get the state of the edge to tree object.
   * An object is either initialized (before commit) or committed (ready to use).
   * \return INITIALIZED or COMMITTED
   */
  inline int
  get_state () const
  {
    return state;
  }

  /* Setter functions */

  /** A single (tree, local edge) value connected to a global edge is added to the edge_to_tree_list.
   * \param [in] cmesh must be committed.
   * \param [in] global_edge_id The global id of the edge to be added.
   * \param [in] ltreeid The local tree id of a tree that \a global_edge_id is connected to.
   * \param [in] tree_edge The local edge id of \a ltreeid that \a global_edge_id is connected to.
   */
  void
  add_edge_to_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_edge_id, const t8_locidx_t ltreeid,
                      const int tree_edge);

  /** Mark as ready for commit. Meaning that all
   * global edge ids have been added.
   * After commit, no edge ids can be added anymore.
   * \param [in] cmesh A committed cmesh to which the global edge ids are associated.
   */
  void
  commit (const t8_cmesh_t cmesh);

  /** Check whether this instance is committed.
    *
    * \return True if committed. Thus all entries have been set.
    */
  inline bool
  is_committed () const
  {
    return state == COMMITTED;
  }

  /** Compare with another instance of this class.
    *
    * \param [in] other   The other list to compare.
    * \return             True if and only if the stored edge indices match.
    */
  inline int
  is_equal (const t8_cmesh_edge_conn_edge_to_tree& other) const
  {
    /* Two instances are equal if and only if their
    * states are equal and the stored edges are equal. */
    return state == other.state && edge_to_tree == other.edge_to_tree;
  }

  /** Equality operator.
    *
    * \param [in] other   The other list to compare.
    * \return             True if and only if the stored edge indices match.
    */
  inline bool
  operator== (const t8_cmesh_edge_conn_edge_to_tree& other) const
  {
    return is_equal (other);
  }

  /** Typedef for the iterator type */
  typedef ett_storage_type::const_iterator const_iterator;

  /** Iterator begin */
  inline const_iterator
  begin () const
  {
    return edge_to_tree.begin ();
  }

  /** Iterator end */
  inline const_iterator
  end () const
  {
    return edge_to_tree.end ();
  }

  friend struct t8_cmesh_edge_connectivity;

 private:
  /** For each global edge id sort the list of
   * (tree_id, tree_edge) pairs according to
   * tree_id and tree_edge index.
   * Example: (1, 3), (0, 0), (1, 0)
   * becomes: (0, 0), (1, 0), (1, 3)
   */
  void
  sort_list_by_tree_id ();

  /** Check that all local trees and edges of a given cmesh are mapped to a global id.
    *
    * \param [in] cmesh    A committed cmesh.
    * \return             True if and only if each local tree and edge of \a cmesh is associated with a global id.
    */
  int
  contains_all_edges (const t8_cmesh_t cmesh) const;

  /** The actual data storage mapping global edge ids to a list
   * local trees and tree edges. */
  ett_storage_type edge_to_tree;

  /** Stores the state of this instance. */
  enum {
    INITIALIZED, /*< Can currently be filled with entries. */
    COMMITTED    /*< Is filled and cannot be changed. */
  } state;
};

#endif /* !T8_CMESH_EDGE_CONN_EDGE_TO_TREE_HXX */
