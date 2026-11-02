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

/** \file t8_forest_splitters.hxx
 * The partition of a forest as seen by one process, in O(log P) memory, and the routing of items to the processes
 * that own their part of the space-filling curve.
 *
 * The ranks [0, P) are bisected recursively: a rank interval [lo, hi) of at least two ranks is split into the lower
 * half [lo, mid) and the upper half [mid, hi) with mid = lo + ceil ((hi - lo) / 2). Each process stores, for the at
 * most ceil (log2 P) intervals that contain it and another rank, where the upper half starts: the global index of the
 * first element of rank mid and the position on the space-filling curve of the first element of [mid, hi). This
 * replaces the element offsets, first descendants and tree offsets of all processes (O(P) each).
 *
 * The splitters are computed in a bottom-up pass over the bisection, in which each process exchanges a summary of its
 * half with a partner in the other half. Items are routed top-down: in each interval, each process sends the items
 * that belong to the other half to its partner. Both passes take ceil (log2 P) rounds with one message per process
 * and round, like the hypercube exchange of the global-sort 2:1 balance in Suh and Isaac, "Evaluation of a minimally
 * synchronous algorithm for 2:1 octree balance" (SC20), but for any number of processes: rank lo + j of the lower
 * half and rank mid + j of the upper half are partners for j < hi - mid. If the lower half has one rank more, its last
 * rank mid - 1 is paired with the last rank hi - 1 of the upper half, which sends it the summary of the upper half and
 * receives its items for the upper half.
 *
 * A position belongs to the last process with elements whose first element is not after it, or to the first process
 * with elements if it lies before all elements. So the processes with elements cover the whole space-filling curve,
 * and processes without elements own nothing.
 */

#ifndef T8_FOREST_SPLITTERS_HXX
#define T8_FOREST_SPLITTERS_HXX

#include <t8.h>
#include <t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <compare>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

/** A position on the space-filling curve of a forest: a descendant at the forest's maxlevel in a tree.
 * Positions are ordered by tree and then by the linear id of the descendant. */
struct t8_forest_position
{
  t8_gloidx_t gtreeid; /**< The global id of the tree. */
  t8_linearidx_t id;   /**< The linear id of the descendant at the maxlevel of the forest. */

  /** Order by tree, then by linear id. */
  friend auto
  operator<=> (const t8_forest_position &, const t8_forest_position &)
    = default;
};

/** A rank interval of the recursive bisection that contains this process, see \ref t8_forest_splitters.hxx */
struct t8_forest_bisection
{
  int lo;                          /**< The first rank of the interval. */
  int mid;                         /**< The first rank of the upper half. */
  int hi;                          /**< One past the last rank of the interval. */
  bool lower_has_elements;         /**< Whether a rank of the lower half [lo, mid) has elements. */
  bool upper_has_elements;         /**< Whether a rank of the upper half [mid, hi) has elements. */
  t8_forest_position mid_position; /**< If \a upper_has_elements, the position of the first element of [mid, hi). */
  t8_gloidx_t mid_element;         /**< The global index of the first element of rank \a mid. */
};

/** The partition of a forest as seen by one process, see \ref t8_forest_splitters.hxx */
struct t8_forest_splitters
{
  int rank;                          /**< The rank of this process. */
  t8_gloidx_t first_element;         /**< The global index of the first local element. */
  t8_gloidx_t global_num_elements;   /**< The number of elements of all processes. */
  bool has_elements;                 /**< Whether this process has elements. */
  t8_forest_position first_position; /**< If \a has_elements, the position of the first local element. */
  t8_forest_position next_position;  /**< The position of the first element of a higher rank, or \a end_position. */
  t8_forest_position end_position;   /**< The position after all trees, (global number of trees, 0). */

  /** The intervals of the bisection that contain this rank and another rank, from [0, P) inwards. */
  std::vector<t8_forest_bisection> levels;
};

/** Compute the splitters of a committed forest and store them in the forest, replacing existing ones.
 * Takes ceil (log2 P) rounds of one message per process.
 * \param [in,out] forest  A committed forest.
 * \note This function is collective.
 */
void
t8_forest_splitters_compute (t8_forest_t forest);

/** Free the splitters of a forest, if it has any.
 * \param [in,out] forest  A forest.
 */
void
t8_forest_splitters_destroy (t8_forest_t forest);

/** The splitters of a committed forest.
 * \param [in] forest  A committed forest.
 * \return The splitters of \a forest.
 */
const t8_forest_splitters &
t8_forest_get_splitters (const t8_forest_t forest);

/** Check the splitters of a forest against its element offsets, first descendants and tree offsets of all processes.
 * The positions are only checked in forests without deleted elements.
 * \param [in] forest  A committed forest with splitters and all three offset arrays.
 * \return True if the splitters match the offset arrays.
 */
int
t8_forest_splitters_match_offsets (const t8_forest_t forest);

/** The positions of the first and last descendant of an element at the maxlevel of the forest. The positions of all
 * descendants of the element lie between them.
 * \param [in] forest   A committed forest.
 * \param [in] gtreeid  The global id of the tree of \a element.
 * \param [in] eclass   The eclass of the tree.
 * \param [in] element  An element of the tree.
 * \param [out] first   The position of the first descendant of \a element.
 * \param [out] last    The position of the last descendant of \a element.
 */
void
t8_forest_element_positions (const t8_forest_t forest, const t8_gloidx_t gtreeid, const t8_eclass_t eclass,
                             const t8_element_t *element, t8_forest_position *first, t8_forest_position *last);

/** The halves of a bisection interval to which an item is routed, see \ref t8_forest_route. */
enum t8_forest_route_halves {
  T8_FOREST_ROUTE_NONE = 0,  /**< The item is dropped. */
  T8_FOREST_ROUTE_LOWER = 1, /**< The item goes to the lower half. */
  T8_FOREST_ROUTE_UPPER = 2, /**< The item goes to the upper half. */
  T8_FOREST_ROUTE_BOTH = 3   /**< A copy of the item goes to each half. */
};

/** The halves of a bisection interval whose processes own a position in [\a first, \a last].
 * \param [in] level   A bisection interval of the splitters of a forest, which owns a position in [first, last].
 * \param [in] first   The first position of the item.
 * \param [in] last    The last position of the item, not before \a first.
 * \param [in] unique  If true, items that are owned by both halves are dropped, so that each item arrives at its
 *                     unique owner, if it has one.
 * \return The halves to route the item to, a \ref t8_forest_route_halves.
 */
int
t8_forest_route_halves_positions (const t8_forest_bisection &level, const t8_forest_position &first,
                                  const t8_forest_position &last, const bool unique);

/** The halves of a bisection interval whose processes own a global element index in [\a first, \a last].
 * \param [in] level   A bisection interval of the splitters of a forest, which owns an index in [first, last].
 * \param [in] first   The first global element index of the item.
 * \param [in] last    The last global element index of the item, not smaller than \a first.
 * \param [in] unique  If true, items that are owned by both halves are dropped.
 * \return The halves to route the item to, a \ref t8_forest_route_halves.
 */
int
t8_forest_route_halves_elements (const t8_forest_bisection &level, const t8_gloidx_t first, const t8_gloidx_t last,
                                 const bool unique);

/** One round of \ref t8_forest_route: send a buffer to the partner of this process in the other half of a bisection
 * interval and receive the buffers of the partners in the other half.
 * \param [in] forest  A committed forest with splitters.
 * \param [in] level   A bisection interval of the splitters of \a forest.
 * \param [in] send    The bytes for the other half.
 * \param [in,out] recv  The received bytes are appended.
 * \note This function is collective over the ranks of \a level.
 */
void
t8_forest_route_exchange (const t8_forest_t forest, const t8_forest_bisection &level, const std::vector<char> &send,
                          std::vector<char> &recv);

/** Route items to the processes that own them: in each bisection interval containing this process, from [0, P)
 * inwards, the items are split by \a halves, and those for the other half are sent to the partner there.
 * \tparam TItem  A trivially copyable type, which is sent as bytes.
 * \param [in] forest     A committed forest with splitters.
 * \param [in,out] items  On input the items of this process, on output the items that were routed to it.
 * \param [in] halves     Called as halves (item, level) with a \ref t8_forest_bisection that owns the item, returns
 *                        the \ref t8_forest_route_halves to route the item to.
 * \param [in] merge      Called as merge (items) after each round, with the kept items followed by the received
 *                        ones. May reorder, combine or remove items, e.g. duplicates.
 * \note This function is collective.
 */
template <typename TItem, typename THalves, typename TMerge>
void
t8_forest_route (const t8_forest_t forest, std::vector<TItem> &items, THalves &&halves, TMerge &&merge)
{
  static_assert (std::is_trivially_copyable_v<TItem>, "Routed items are sent as bytes");
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  std::vector<TItem> keep;
  std::vector<char> send, recv;
  for (const t8_forest_bisection &level : splitters.levels) {
    const int own_half = splitters.rank < level.mid ? T8_FOREST_ROUTE_LOWER : T8_FOREST_ROUTE_UPPER;
    keep.clear ();
    send.clear ();
    for (const TItem &item : items) {
      const int item_halves = halves (item, level);
      if (item_halves & own_half) {
        keep.push_back (item);
      }
      if (item_halves & ~own_half & T8_FOREST_ROUTE_BOTH) {
        const char *bytes = reinterpret_cast<const char *> (&item);
        send.insert (send.end (), bytes, bytes + sizeof (TItem));
      }
    }
    recv.clear ();
    t8_forest_route_exchange (forest, level, send, recv);
    T8_ASSERT (recv.size () % sizeof (TItem) == 0);
    const size_t num_kept = keep.size ();
    keep.resize (num_kept + recv.size () / sizeof (TItem));
    if (!recv.empty ()) {
      std::memcpy (static_cast<void *> (keep.data () + num_kept), recv.data (), recv.size ());
    }
    items.swap (keep);
    merge (items);
  }
}

/** Route items to the processes that own them, see \ref t8_forest_route.
 * \param [in] forest     A committed forest with splitters.
 * \param [in,out] items  On input the items of this process, on output the items that were routed to it.
 * \param [in] halves     Called as halves (item, level), returns the \ref t8_forest_route_halves of the item.
 * \note This function is collective.
 */
template <typename TItem, typename THalves>
void
t8_forest_route (const t8_forest_t forest, std::vector<TItem> &items, THalves &&halves)
{
  t8_forest_route (forest, items, std::forward<THalves> (halves), [] (std::vector<TItem> &) {});
}

#endif /* !T8_FOREST_SPLITTERS_HXX */
