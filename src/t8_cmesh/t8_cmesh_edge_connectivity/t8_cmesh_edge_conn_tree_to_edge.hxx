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

/** \file t8_cmesh_edge_conn_tree_to_edge.hxx
 * Class to save data structure for tree_to_edge_lists of the cmesh.
 * When the cmesh stores global edge numbers, we require a lookup that
 * matches a tree and its local edge to a global edge id.
 * This lookup is encoded in the t8_cmesh_edge_conn_tree_to_edge struct.
 */

/* TODO:
 *  It is probably best to set all global ids of a single tree as one attribute.
 * That way we can store it as a single arrays of id's.
 * We will probably most often want to access all ids of one tree, so a function returning an array of ids
 * is a good idea anyway and if we already store them as such then we do not need to do any data movement
 * when accessing.
 *
 * On the downside we will only have a "set all ids of a tree" function and no "set this single id for this tree and this edge" function.
 */

#ifndef T8_CMESH_EDGE_CONN_TREE_TO_EDGE_HXX
#define T8_CMESH_EDGE_CONN_TREE_TO_EDGE_HXX

#include <algorithm>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_conn_edge_to_tree.hxx>
#include <span>

/** forward declaration of the class needed since the two class headers include each other. */
struct t8_cmesh_edge_conn_edge_to_tree;

/**
 * A class to hold the tree to edge connectivity of a cmesh.
 */
struct t8_cmesh_edge_conn_tree_to_edge
{
 public:
  /** Standard constructor. Does nothing. */
  t8_cmesh_edge_conn_tree_to_edge (): current_state (state::EMPTY)
  {
  }

  /** Constructor from a cmesh where all the attributes are set.
   * Currently unclear if we implement this eventually.
   * If we do so: Should the cmesh be already committed, or in pre-commit state but attributes set?
   *
   * \note This function is not implemented yet.
   */
  t8_cmesh_edge_conn_tree_to_edge ([[maybe_unused]] const t8_cmesh_t cmesh)
  {
    // TODO: Remove the [[maybe unused]] qualifier when implemented
    SC_ABORT ("not implemented.");
  }

  /** Constructor from a cmesh and a given edge to tree connectivity.
   *
   * \param [in] cmesh_from A committed cmesh.
   * \param [in] cmesh      An initialized but not committed cmesh that is to be derived from \a cmesh_from.
   * \param [in] ett        A committed edge to tree connectivity for \a cmesh_from.
   *
   * As a result a tree to edge connectivity for \a cmesh will be constructed.
   * \note \a cmesh_from must be committed.
   * \note \a cmesh must not be committed.
   * \note \a ett must be committed.
   * \note This does not work until issue #923 https://github.com/DLR-AMR/t8code/issues/923 is resolved.
   */
  t8_cmesh_edge_conn_tree_to_edge (const t8_cmesh_t cmesh_from, const t8_cmesh_t cmesh,
                                       const struct t8_cmesh_edge_conn_edge_to_tree &ett);

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
    T8_ASSERT (num_edges >= 0);
    T8_ASSERT (global_tree_edges != NULL);

    /* TODO: we currently do not check whether the num_edges argument
    *       matches the number of edges of the tree.
    *       We cannot do it here, since this function call happens before commit,
    *       thus we might not even know the eclass of the tree.
    *       Maybe it is possible to check this during t8_cmesh_add_attributes?
    */

    /* We copy the data directly, hence set data_persiss to 0 */
    const int data_persists = 0;
    t8_debugf ("Setting %i global edges for global tree %" T8_GLOIDX_FORMAT ".\n", num_edges, global_tree);
    t8_cmesh_set_attribute_gloidx_array (cmesh, global_tree, t8_get_package_id (),
                                         T8_CMESH_GLOBAL_EDGES_ATTRIBUTE_KEY, global_tree_edges, num_edges,
                                         data_persists);
    current_state = state::FILLED;
  }

  /* TODO: What if the attribute is not set? error handling */
  /** Return the global edge indices of a local tree.
   * \param [in] cmesh A committed cmesh.
   * \param [in] local_tree A local tree in \a cmesh.
   * \return An array of length \a num_edges containing the global edge ids of \a local_tree's edges.
  */
  inline const std::span<const t8_gloidx_t>
  get_global_edges (const t8_cmesh_t cmesh, const t8_locidx_t local_tree) const
  {
    T8_ASSERT (t8_cmesh_is_committed (cmesh));

    /* Get num tree_edges to create a view */
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, local_tree);
    const int num_tree_edges = t8_eclass_num_edges[tree_class];

    const t8_gloidx_t *global_edges = t8_cmesh_get_attribute_gloidx_array (
      cmesh, t8_get_package_id (), T8_CMESH_GLOBAL_EDGES_ATTRIBUTE_KEY, local_tree, num_tree_edges);
    T8_ASSERT (global_edges != NULL);
    const std::span<const t8_gloidx_t> view (global_edges, num_tree_edges);
    return view;
  }

  /* TODO: What if the attribute is not set? error handling */
  /** Return a single global edge id of a single local edge.
   *
   *
   * \param [in] cmesh A committed cmesh.
   * \param [in] local_tree A local tree of \a cmesh.
   * \param [in] local_tree_edge A local edge of \a local_tree
   * \return The global id of the local edge \a local_tree_edge of \a local_tree.
   */
  t8_gloidx_t
  get_global_edge (const t8_cmesh_t cmesh, const t8_locidx_t local_tree, const int local_tree_edge) const
  {
    T8_ASSERT (t8_cmesh_is_committed (cmesh));

    /* Verify that local_tree_edge is in fact a local edge of the tree */
    /* Note: We only perform this check in debugging mode.
    *       In non-debugging mode, using a edge index beyond the trees index allows
    *       for a potential attacker to gain access to memory possibly not owned by the caller.
    *       We do not check in non-debugging mode for (obvious) performance reasons. */
    T8_ASSERT (0 <= local_tree_edge);
    const std::span<const t8_gloidx_t> global_edges = get_global_edges (cmesh, local_tree);
    T8_ASSERT (global_edges.size () > static_cast<size_t> (local_tree_edge));
    return global_edges[local_tree_edge];
  }

  friend struct t8_cmesh_edge_connectivity;

 private:
  enum struct state {
    EMPTY, /*< Is initialized but empty. */
    FILLED /*< Is filled with at least one entry. */
  };

  /** Return the state of this object. */
  inline state
  get_state () const
  {
    return current_state;
  }

  state current_state;
};

#endif /* !T8_CMESH_EDGE_CONN_TREE_TO_EDGE_HXX */
