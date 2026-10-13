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

#ifndef T8_CMESH_REINDEX_TREES_H
#define T8_CMESH_REINDEX_TREES_H

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <map>

/**
 * Compute a spatially coherent reindexing of the coarse mesh trees.
 *
 * Constructs an auxiliary bounding box forest containing the centers of the
 * original coarse mesh trees. The forest is adaptively refined until each
 * leaf element contains at most one tree center. The leaves are then traversed
 * in SFC order to assign new indices to the original trees.
 *
 * \param [in] cmesh  The original coarse mesh whose trees are to be reindexed.
 * \param [in] comm   MPI communicator used to construct the auxiliary forest.
 *
 * \return A map from the original global tree IDs to their new SFC-based indices.
 *
 * \note This function only computes the reindexing and does not modify the
 *       original coarse mesh.
 */

std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm = sc_MPI_COMM_SELF);

/**
 * Apply a tree reindexing to the coarse mesh stash in place.
 *
 * Updates the tree IDs stored in the stash according to the given mapping,
 * allowing the coarse mesh to use the newly computed tree ordering.
 *
 * \param [in,out] stash         The coarse mesh stash to be modified.
 * \param [in]     tree_reindex  Mapping from original global tree IDs to
 *                              their new indices.
 */
void
t8_cmesh_tree_perform_reindex_inplace (t8_stash_t &stash, const std::map<t8_gloidx_t, t8_gloidx_t> &tree_reindex);

#endif /* !T8_CMESH_REINDEX_TREES_H */
