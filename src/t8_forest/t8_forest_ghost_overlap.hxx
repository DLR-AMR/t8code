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

/** \file t8_forest_ghost_overlap.hxx
 * A ghost layer of a given width: the leaves of other ranks in a region around the part of this rank, for
 * overlapping domain decomposition (\ref t8_forest_set_ghost_overlap).
 *
 * The region of a rank is defined at a level l and a number of layers m: the elements of level l whose distance from
 * the part of the rank, counted in steps between level-l elements with a common vertex (inside the trees and across
 * tree boundaries), is at most m. Distance 0 are the level-l elements intersecting the part. So the width of the
 * overlap is at least m times the size of the level-l elements, independent of the refinement of the leaves. The
 * region is built from the cover elements of the part (the maximal elements between its first and its next position,
 * \ref t8_forest_touching_ranks.hxx): a finer one is replaced by its ancestor of level l, a coarser one by its
 * descendants of level l that touch its boundary (the other descendants and their neighbors lie in the part). Then the
 * region grows m times by the elements of level l that share a vertex with an element of the region. With m >= 1, the
 * region contains every leaf of another rank that shares a vertex with a local leaf, so the ghost layer contains the
 * vertex ghost layer.
 *
 * The overlap ghosts of a rank are the leaves of other ranks that intersect an element of its region: the leaves
 * inside a region element and the leaf containing one. The region elements that do not lie in the own part are
 * routed to the ranks that own a part of them (\ref t8_forest_route, O(log P) rounds), which collect their leaves
 * intersecting them. The relation is symmetric: if a leaf of rank B intersects a level-l element X at distance at most
 * m from an element Y intersecting the part of rank A, then X intersects the part of B and Y lies in the region of B.
 * The ghost exchange relies on this (it receives from the ranks it sends to). To be safe, every rank that sends leaves
 * to another one also routes a notice to the owner of its first position, and the senders of overlap ghosts are
 * remote ranks of the ghost layer even if they receive nothing (an empty message).
 */

#ifndef T8_FOREST_GHOST_OVERLAP_HXX
#define T8_FOREST_GHOST_OVERLAP_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_touching_ranks.hxx>
#include <vector>

/** The local leaves in the overlap regions of other ranks, and the ranks that have leaves in the region of this rank,
 * see \ref t8_forest_ghost_overlap.hxx.
 * \param [in,out] forest   A committed forest with splitters. The boundary connectivity of its coarse mesh is built if
 *                          necessary.
 * \param [in]     level    The level l of the region elements.
 * \param [in]     layers   The number of layers m around the part of this rank.
 * \param [out]    mirrors  The local leaves with the ranks whose region they intersect, sorted by rank, tree and
 *                          index, without duplicates.
 * \param [out]    senders  The ranks other than this one that own a leaf intersecting the region of this rank, sorted.
 * \param [out]    ranges   The ranges of the ranks other than this one that send overlap ghosts to this rank or
 *                          receive them from it, sorted by rank. With them, the owners of the overlap ghosts are
 *                          known (\ref t8_forest_touching_ranks_add).
 * \note This function is collective.
 */
void
t8_forest_ghost_overlap_mirrors (t8_forest_t forest, int level, int layers, std::vector<t8_forest_mirror> &mirrors,
                                 std::vector<int> &senders, std::vector<t8_forest_rank_range> &ranges);

#endif /* !T8_FOREST_GHOST_OVERLAP_HXX */
