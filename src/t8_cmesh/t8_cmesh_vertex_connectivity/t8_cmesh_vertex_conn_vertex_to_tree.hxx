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

/** \file t8_cmesh_vertex_conn_vertex_to_tree.hxx
 * Class to save data structures cmesh for vertex_to_tree_lists
 */

#ifndef T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX
#define T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX

#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity_types.hxx>

/** forward declaration of ttv class needed since the two class headers include each other. */
class t8_cmesh_vertex_conn_tree_to_vertex;

/** This class stores the vertex to tree lookup for
 * global vertex indices for a cmesh.
 * Thus, given a global vertex id the class provides
 * information about the trees the vertex belongs to and
 * the corresponding local vertex ids inside these trees.
 *
 * In particular, this class stores the lookup
 *
 * global_vertex_id -> List of (tree, tree_local_vertex)
 *
 * for a cmesh.
 * It is the opposite lookup as \ref t8_cmesh_vertex_conn_tree_to_vertex.hxx
 *
 * The global vertex ids must not be contiguous, that is, we have some set
 *
 * {I_0 < I_1 < ...< I_N} of natural numbers corresponding to the N+1 vertices.
 *
 * I_0 does not have to be 0 and I_N does not have to be N.
 *
 * So we need lookup: I_i -> i
 *  store this in a hash table.
 *
 * datatypes:
 *
 * global id: t8_gloidx_t
 * (tree_id, tree_vertex): std::pair<t8_locidx_t, int> = TV_PAIR
 * List of (tree_id, tree_vertex): std::vector<PAIR> = TV_LIST
 * Table global_id -> TV_LIST: std::unordered_map<t8_gloidx_t, TV_LIST>
 *
 */
class t8_cmesh_vertex_conn_vertex_to_tree {
 public:
  /** Standard constructor.
   * Initializes the class and allows setting vertex entries
   * via \ref add_vertex_to_tree
   */
  t8_cmesh_vertex_conn_vertex_to_tree (): state (INITIALIZED)
  {
  }

  /** Function to fill vtt from cmesh and ttv information.
   * Sets all global ids and associated tree vertices from
   * the given input class.
   * Afterwards, the class is set to committed and can be used.
   *
   * \param [in] cmesh A committed cmesh with set tree to vertex entries.
   * \param [in] ttv A filled tree to vertex list for \a cmesh.
  */
  void
  build_from_ttv (const t8_cmesh_t cmesh, t8_cmesh_vertex_conn_tree_to_vertex& ttv);

  /* Setter functions */
  /** Given a cmesh, build up the vertex_to_tree.
   * \param [in] cmesh An initialized but not yet committed cmesh.
   * The cmesh must not be committed, but all tree information and neighbor information must
   * have been set.
   * Currently, \a cmesh has to be replicated. */
  void
  set_vertex_to_tree_list (const t8_cmesh_t cmesh);

  /** Get the list of global trees and local vertex ids a global vertex is connected to.
   *
   * \param [in] global_vertex_id The global id of a vertex in the cmesh.
   * \return The list of local tree ids and local vertex ids of \a global_vertex_id.
   */
  inline const tree_vertex_list&
  get_tree_list_of_vertex (const t8_gloidx_t global_vertex_id) const
  {
    T8_ASSERT (is_committed ());
    T8_ASSERT (0 <= global_vertex_id);

    /* Use at() to look for the vertex entry.
    * If the entry does not exist an exception of
    * type std::out_of_range is thrown. */

    try {
      return vertex_to_tree.at (global_vertex_id);
    } catch (const std::out_of_range& e) {
      t8_errorf ("ERROR: Could not find vertex %li for cmesh.\n", static_cast<long>(global_vertex_id));
      SC_ABORTF ("Caught exception 'out of range': %s\n", e.what ());
    }
  }

  /** Get the state of the vertex to tree object.
   * An object is either initialized (before commit) or committed (ready to use).
   * \return INITIALIZED or COMMITTED
   */
  inline int
  get_state () const
  {
    return state;
  }

  /* Setter functions */

  /** A single (tree, local vertex) value connected to a global vertex is added to the vertex_to_tree_list.
   * \param [in] cmesh must be committed.
   * \param [in] global_vertex_id The global id of the vertex to be added.
   * \param [in] ltreeid The local tree id of a tree that \a global_vertex_id is connected to.
   * \param [in] tree_vertex The local vertex id of \a ltreeid that \a global_vertex_id is connected to.
   */
  void
  add_vertex_to_tree (const t8_cmesh_t cmesh, const t8_gloidx_t global_vertex_id, const t8_locidx_t ltreeid,
                      const int tree_vertex);

  /** Mark as ready for commit. Meaning that all
   * global vertex ids have been added.
   * After commit, no vertex ids can be added anymore.
   * \param [in] cmesh A committed cmesh to which the global vertex ids are associated.
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
    * \return             True if and only if the stored vertex indices match.
    */
  inline int
  is_equal (const t8_cmesh_vertex_conn_vertex_to_tree& other) const
  {
    /* Two instances are equal if and only if their
    * states are equal and the stored vertices are equal. */
    return state == other.state && vertex_to_tree == other.vertex_to_tree;
  }

  /** Equality operator.
    *
    * \param [in] other   The other list to compare.
    * \return             True if and only if the stored vertex indices match.
    */
  inline bool
  operator== (const t8_cmesh_vertex_conn_vertex_to_tree& other) const
  {
    return is_equal (other);
  }

  /** Typedef for the iterator type */
  typedef vtt_storage_type::const_iterator const_iterator;

  /** Iterator begin */
  inline const_iterator
  begin () const
  {
    return vertex_to_tree.begin ();
  }

  /** Iterator end */
  inline const_iterator
  end () const
  {
    return vertex_to_tree.end ();
  }

  friend struct t8_cmesh_vertex_connectivity;

 private:
  /** For each global vertex id sort the list of
   * (tree_id, tree_vertex) pairs according to
   * tree_id and tree_vertex index.
   * Example: (1, 3), (0, 0), (1, 0)
   * becomes: (0, 0), (1, 0), (1, 3)
   */
  void
  sort_list_by_tree_id ();

  /** Check that all local trees and vertices of a given cmesh are mapped to a global id.
    *
    * \param [in] cmesh    A committed cmesh.
    * \return             True if and only if each local tree and vertex of \a cmesh is associated with a global id.
    */
  int
  contains_all_vertices (const t8_cmesh_t cmesh) const;

  /** The actual data storage mapping global vertex ids to a list
   * local trees and tree vertices. */
  vtt_storage_type vertex_to_tree;

  /** Stores the state of this instance. */
  enum {
    INITIALIZED, /*< Can currently be filled with entries. */
    COMMITTED    /*< Is filled and cannot be changed. */
  } state;
};

#endif /* !T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX */
