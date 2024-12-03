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

/** \file t8_cmesh_vertex_conn_vertex_to_tree.hxx
 * Class to save data structure for vertex_to_tree_lists
 */

#ifndef T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX
#define T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX

#include <vector>
#include <unordered_map>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_conn_tree_to_vertex.hxx>

/* forward declaration of ttv class needed since the two class headers include each other. */
class t8_cmesh_vertex_conn_tree_to_vertex;

/*
 *  notes during development
 * 
 * This class stores the lookup
 * 
 * global_vertex_id -> List of (tree, tree_local_vertex) 
 * 
 * for a cmesh.
 * It is the opposite lookup as t8_cmesh_vertex_conn_tree_to_vertex
 * 
 * The global vertex ids must not be contiguous, that is, we have some set
 * 
 * {I_0 < I_1 < ...< I_N} of natural numbers corresponding to the N+1 vertices.
 * 
 * I_0 does not have to be 0 and I_N does not have to be N.
 * 
 * 
 * So we need lookup: I_i -> i
 *  store this in a hash table.
 * 
 * dattypes:
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

  /* Variable type for (tree_id, tree_vertex_id) pair */
  using tree_vertex_pair = std::pair<t8_locidx_t, int>;

  /* list of tree vertex pairs, each global vertex id maps to 
    * such a list. */
  using tree_vertex_list = std::vector<tree_vertex_pair>;

  using vtt_storage_type = std::unordered_map<t8_gloidx_t, tree_vertex_list>;

  /* Setter functions */
  /* Given a cmesh, build up the vertex_to_tree.
   * \return: some error value to be specified.
   * The cmesh must not be committed, but all tree information and neighbor information must
   * have been set. 
   * Currently, \a cmesh has to be replicated. */
  void
  set_vertex_to_tree_list (const t8_cmesh_t cmesh);

  const tree_vertex_list&
  get_tree_list_of_vertex (t8_gloidx_t global_vertex_id) const;

  const int
  get_state ()
  {
    return state;
  }

  /* Setter functions */
  /* A single value is added to the vertex_to_tree_list.
   * \a cmesh must be committed. */
  void
  add_vertex_to_tree (const t8_cmesh_t cmesh, t8_gloidx_t global_vertex_id, t8_locidx_t ltreeid, int tree_vertex);

  /* Mark as ready for commit. Meaning that all 
   * global vertex ids have been added.
   * After commit, no vertex ids can be added anymore. */
  void
  commit (const t8_cmesh_t cmesh);

  /**
   * @brief Check whether this instance is committed.
   * 
   * @return int True if committed. Thus all entries have been set.
   */
  int
  is_committed () const;

  /**
   * @brief Compare with another instance of this class.
   * 
   * @param other 
   * @return int True if and only if the stored vertex indices match.
   */
  int
  is_equal (const t8_cmesh_vertex_conn_vertex_to_tree& other) const;

  /**
   * @brief Equality operator. Implement
   * 
   * @param other 
   * @return true 
   * @return false 
   */
  bool
  operator== (const t8_cmesh_vertex_conn_vertex_to_tree& other) const;

  /** Typedef for the iterator type */
  typedef vtt_storage_type::const_iterator const_iterator;

  const_iterator
  begin () const
  {
    return vertex_to_tree.begin ();
  }
  const_iterator
  end () const
  {
    return vertex_to_tree.end ();
  }

  friend struct t8_cmesh_vertex_connectivity;

 private:
  /* For each global vertex id sort the list of
   * (tree_id, tree_vertex) pairs according to
   * tree_id and tree_vertex index.
   * Example: (1, 3), (0, 0), (1, 0)
   * becomes: (0, 0), (1, 0), (1, 3)
   */
  void
  sort_list_by_tree_id ();

  /**
   * @brief Check that all local trees and vertices of a given cmesh are mapped to a global id.
   * 
   * @param cmesh A committed cmesh.
   * @return int True if and only if each local tree and vertex of \a cmesh is associated with a global id.
   */
  int
  contains_all_vertices (const t8_cmesh_t cmesh) const;

  /* The actual data storage mapping global vertex ids to a list
   * local trees and tree vertices. */
  vtt_storage_type vertex_to_tree;

  /** Stores the state of this instance. */
  enum {
    INITIALIZED, /*< Can currently be filled with entries. */
    COMMITTED    /*< Is filled and cannot be changed. */
  } state;
};

#endif /* !T8_CMESH_VERTEX_CONN_VERTEX_TO_TREE_HXX */
