/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_tree_reindex.hxx
 * Internal declarations for reindexing coarse trees by geometric locality.
 */

#ifndef T8_CMESH_TREE_REINDEX_HXX
#define T8_CMESH_TREE_REINDEX_HXX

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

#include <map>

/** Calculate the old-to-new coarse-tree ID mapping from the auxiliary forest's SFC order.
 *
 * The tree centers are computed from the vertex attributes of all coarse trees.
 * The auxiliary forest is adapted until its leaves isolate the centers or reach
 * the maximum refinement level. Multiple centers in one terminal leaf are ordered
 * by their original global tree IDs.
 *
 * \param [in] cmesh  Uncommitted coarse mesh with a complete replicated stash.
 * \param [in] comm   Communicator used for the auxiliary forest. Defaults to MPI self.
 * \return            Map from original global tree ID to its new SFC-based ID.
 *
 * \note This function only calculates the permutation and does not modify \a cmesh.
 */
std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_calculate_tree_reindexing_map (t8_cmesh_t cmesh, sc_MPI_Comm comm = sc_MPI_COMM_SELF);

/** Apply an already calculated old-to-new tree-ID permutation to a stash.
 * Changes tree classes, attributes, and both endpoints of every face join, then
 * restores the canonical ordering of the stash entries.
 *
 * \param [in,out] stash         Stash containing the tree records to reindex.
 * \param [in]     tree_reindex  Map from old global IDs to new global IDs.
 */
void
t8_cmesh_tree_perform_reindex_inplace (t8_stash_t &stash, const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex);

/** Calculate and apply SFC-based coarse-tree reindexing.
 *
 * Calls \ref t8_cmesh_calculate_tree_reindexing_map and
 * \ref t8_cmesh_tree_perform_reindex_inplace. If requested with
 * \ref t8_cmesh_enable_tree_reordering, the original-to-new mapping is
 * retained in the coarse mesh for retrieval after commit.
 *
 * \param [in,out] cmesh  Uncommitted coarse mesh with a complete replicated stash.
 * \note Currently not supported for partitioned coarse meshes.
 */
void
t8_cmesh_reindex_trees (t8_cmesh_t cmesh);

#endif /* T8_CMESH_TREE_REINDEX_HXX */
