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

/** \file t8_forest_touching_ranks.hxx
 * The ranks whose part of a forest touches the part of this rank, together with the positions they own.
 *
 * The part of a rank with elements is the union of the elements between its first position and the first position of
 * the next rank with elements (\ref t8_forest_splitters.hxx). It is covered by the maximal elements in this range,
 * at most 2^d - 1 per level at both ends of the range plus the trees in between. For each of these cover elements and
 * each of its boundaries of dimension at least min_dim, the neighbors of the same level that do not lie in the own
 * range are routed to the ranks that own a part of them (\ref t8_forest_route), together with the boundary of the
 * neighbor and the own range. A receiving rank checks whether its part touches this boundary with dimension at least
 * min_dim, recursively over the children touching the boundary. If so, the sender is a touching rank. Touching is
 * symmetric, so every rank learns all ranks that touch it in one routing pass.
 *
 * The parts of the ranks do not change when a forest is adapted or balanced without repartitioning (unless elements
 * are removed), so such forests inherit the touching ranks of the forest they are derived from.
 *
 * The owner of a position in the range of this rank or of a touching rank is known. Between these ranges, only the
 * ranks in between are known as bounds for the owner. Every rank that owns a leaf touching a local leaf with dimension
 * at least min_dim touches this rank, so the owner searches of the ghost layer only need owners in known ranges.
 *
 * The same neighbors give the ghost layer (\ref t8_forest_touching_ranks_mirrors): a leaf of another rank touches the
 * part of the sender if and only if it contains a neighbor or lies in one and touches its boundary shared with the
 * cover element with dimension at least min_dim. So the receivers collect these leaves, which takes time proportional
 * to their number times the depth, instead of searching the owners at all boundaries of all leaves near the
 * boundary of the own part.
 */

#ifndef T8_FOREST_TOUCHING_RANKS_HXX
#define T8_FOREST_TOUCHING_RANKS_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_splitters.hxx>
#include <vector>

/** The positions owned by a rank: from its first element to the first element of the next rank with elements. */
struct t8_forest_rank_range
{
  int rank;                 /**< The rank. */
  t8_forest_position first; /**< The position of the first element of the rank. */
  t8_forest_position next;  /**< The position of the first element of the next rank with elements, or the end. */
};

/** The ranks whose part of a forest touches the part of this rank, see \ref t8_forest_touching_ranks.hxx */
struct t8_forest_touching_ranks
{
  int min_dim; /**< The parts touch with dimension at least min_dim. */

  /** The ranges of this rank (if it has elements) and of all touching ranks, sorted by rank. */
  std::vector<t8_forest_rank_range> ranges;
};

/** A local leaf that touches the part of another rank, so it is a ghost of that rank. */
struct t8_forest_mirror
{
  int rank;            /**< The rank whose part the leaf touches. */
  t8_locidx_t ltreeid; /**< The local tree of the leaf. */
  t8_locidx_t index;   /**< The index of the leaf in its tree. */

  /** Order by rank, tree and index. */
  friend auto
  operator<=> (const t8_forest_mirror &, const t8_forest_mirror &)
    = default;
};

/** Compute the ranks whose part of a forest touches the part of this rank with dimension at least \a min_dim and
 * store them in the forest, replacing existing ones. Takes one \ref t8_forest_route.
 * \param [in,out] forest   A committed forest.
 * \param [in]     min_dim  The minimal dimension of the contact. If it is the dimension of the forest minus one, only
 *                          faces are used, otherwise the scheme must implement boundaries (standalone).
 * \note This function is collective.
 */
void
t8_forest_touching_ranks_compute (t8_forest_t forest, const int min_dim);

/** Make sure that a forest has the ranks touching it with dimension at least \a min_dim, computing them if it has
 * none or only those touching with a larger dimension.
 * \param [in,out] forest   A committed forest.
 * \param [in]     min_dim  The minimal dimension of the contact.
 * \note This function is collective.
 */
void
t8_forest_touching_ranks_require (t8_forest_t forest, const int min_dim);

/** The local leaves that touch the parts of other ranks with dimension at least \a min_dim, which are the remote
 * elements of the ghost layer. As in \ref t8_forest_touching_ranks_compute, every rank sends the neighbors of its
 * cover elements to the ranks that own a part of them, and these collect their leaves touching the boundary of the
 * neighbor. If the forest has touching ranks for \a min_dim, the neighbors are sent directly to the touching ranks
 * that own a part of them in one round. Otherwise they are routed, and the touching ranks are computed as well.
 * \param [in,out] forest   A committed forest.
 * \param [in]     min_dim  The minimal dimension of the contact, see \ref t8_forest_touching_ranks_compute.
 * \param [out]    mirrors  The leaves with the ranks they touch, sorted by rank, tree and index, without duplicates.
 * \note This function is collective.
 */
void
t8_forest_touching_ranks_mirrors (t8_forest_t forest, const int min_dim, std::vector<t8_forest_mirror> &mirrors);

/** Copy the touching ranks of a forest to a forest with the same partition.
 * \param [in,out] forest  A forest whose ranks own the same positions as those of \a from.
 * \param [in]     from    A forest, which may have touching ranks.
 */
void
t8_forest_touching_ranks_copy (t8_forest_t forest, const t8_forest_t from);

/** Free the touching ranks of a forest, if it has any.
 * \param [in,out] forest  A forest.
 */
void
t8_forest_touching_ranks_destroy (t8_forest_t forest);

/** Check that the range of this rank in the touching ranks of a forest is its range in the splitters.
 * \param [in] forest  A committed forest with splitters and touching ranks.
 * \return True if the ranges match.
 */
int
t8_forest_touching_ranks_match_splitters (const t8_forest_t forest);

/** Bounds for the owner of a position: the owner if the position lies in the range of this rank or of a touching
 * rank, otherwise the ranks between the touching ranks before and after the position.
 * \param [in] forest     A committed forest with touching ranks.
 * \param [in] position   A position.
 * \param [out] lower     A lower bound for the owner of \a position.
 * \param [out] upper     An upper bound for the owner of \a position, equal to \a lower if the owner is known.
 */
void
t8_forest_touching_ranks_owner_bounds (const t8_forest_t forest, const t8_forest_position &position, int *lower,
                                       int *upper);

#endif /* !T8_FOREST_TOUCHING_RANKS_HXX */
