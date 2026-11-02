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

#include <t8_forest/t8_forest_splitters.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_data/t8_shmem.h>
#include <t8_schemes/t8_scheme.hxx>
#include <climits>

/** The summary of a rank interval that is exchanged in the bottom-up pass of \ref t8_forest_splitters_compute. */
struct t8_forest_splitters_summary
{
  t8_gloidx_t num_elements;          /**< The number of elements of the ranks of the interval. */
  t8_forest_position first_position; /**< If num_elements > 0, the position of the first element of the interval. */
};

void
t8_forest_splitters_compute (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_forest_splitters_destroy (forest);
  t8_forest_splitters *splitters = forest->splitters = new t8_forest_splitters;
  const int rank = forest->mpirank;
  splitters->rank = rank;
  splitters->end_position = { forest->global_num_trees, 0 };

  /* The intervals that contain this rank and another, from [0, P) inwards */
  for (int lo = 0, hi = forest->mpisize; hi - lo > 1;) {
    const int mid = lo + (hi - lo + 1) / 2;
    splitters->levels.push_back ({ lo, mid, hi, false, false, {}, 0 });
    (rank < mid ? hi : lo) = mid;
  }

  /* The summary of this rank */
  t8_forest_splitters_summary summary = { forest->local_num_leaf_elements, {} };
  splitters->has_elements = summary.num_elements > 0;
  if (splitters->has_elements) {
    /* The first descendant of the first leaf of the first tree with leaves */
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    t8_locidx_t itree = 0;
    while (t8_forest_get_tree_num_leaf_elements (forest, itree) == 0) {
      itree++;
    }
    const t8_tree_t tree = t8_forest_get_tree (forest, itree);
    summary.first_position = { forest->first_local_tree + itree,
                               scheme->element_get_linear_id (tree->eclass, tree->first_desc, forest->maxlevel) };
    splitters->first_position = summary.first_position;
  }

  /* Bottom-up: exchange the summaries of the halves of each interval and store where the upper half starts */
  std::vector<t8_gloidx_t> lower_num_elements (splitters->levels.size ());
  for (size_t ilevel = splitters->levels.size (); ilevel-- > 0;) {
    t8_forest_bisection &level = splitters->levels[ilevel];
    const bool is_lower = rank < level.mid;
    const int j = rank - (is_lower ? level.lo : level.mid);
    const int num_upper = level.hi - level.mid;
    int targets[2], num_targets = 0, source;
    if (is_lower) {
      /* The last rank of a lower half with one rank more gets the summary of the upper half from its last rank */
      source = j < num_upper ? level.mid + j : level.hi - 1;
      if (j < num_upper) {
        targets[num_targets++] = source;
      }
    }
    else {
      source = targets[num_targets++] = level.lo + j;
      if (j == num_upper - 1 && level.mid - level.lo > num_upper) {
        targets[num_targets++] = level.mid - 1;
      }
    }
    sc_MPI_Request requests[2];
    for (int itarget = 0; itarget < num_targets; itarget++) {
      const int mpiret = sc_MPI_Isend (&summary, sizeof (summary), sc_MPI_BYTE, targets[itarget],
                                       T8_MPI_FOREST_SPLITTERS, forest->mpicomm, requests + itarget);
      SC_CHECK_MPI (mpiret);
    }
    t8_forest_splitters_summary other;
    int mpiret = sc_MPI_Recv (&other, sizeof (other), sc_MPI_BYTE, source, T8_MPI_FOREST_SPLITTERS, forest->mpicomm,
                              sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Waitall (num_targets, requests, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    const t8_forest_splitters_summary lower = is_lower ? summary : other;
    const t8_forest_splitters_summary upper = is_lower ? other : summary;
    level.lower_has_elements = lower.num_elements > 0;
    level.upper_has_elements = upper.num_elements > 0;
    level.mid_position = upper.first_position;
    lower_num_elements[ilevel] = lower.num_elements;
    summary.num_elements = lower.num_elements + upper.num_elements;
    summary.first_position = level.lower_has_elements ? lower.first_position : upper.first_position;
  }
  splitters->global_num_elements = summary.num_elements;
  T8_ASSERT (splitters->global_num_elements == forest->global_num_leaf_elements);

  /* Top-down: the global index of the first element of each upper half */
  t8_gloidx_t first_element = 0;
  for (size_t ilevel = 0; ilevel < splitters->levels.size (); ilevel++) {
    t8_forest_bisection &level = splitters->levels[ilevel];
    level.mid_element = first_element + lower_num_elements[ilevel];
    if (rank >= level.mid) {
      first_element = level.mid_element;
    }
  }
  splitters->first_element = first_element;

  /* The ranks after this one are the upper halves of the intervals in whose lower half this rank is, from the
   * innermost interval outwards. The first of them with elements has the first element after this rank. */
  splitters->next_position = splitters->end_position;
  for (size_t ilevel = splitters->levels.size (); ilevel-- > 0;) {
    const t8_forest_bisection &level = splitters->levels[ilevel];
    if (rank < level.mid && level.upper_has_elements) {
      splitters->next_position = level.mid_position;
      break;
    }
  }
}

void
t8_forest_splitters_destroy (t8_forest_t forest)
{
  delete forest->splitters;
  forest->splitters = NULL;
}

const t8_forest_splitters &
t8_forest_get_splitters (const t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (forest->splitters != NULL);
  return *forest->splitters;
}

int
t8_forest_splitters_match_offsets (const t8_forest_t forest)
{
  T8_ASSERT (forest->splitters != NULL);
  T8_ASSERT (forest->element_offsets != NULL && forest->global_first_desc != NULL && forest->tree_offsets != NULL);
  const t8_forest_splitters &splitters = *forest->splitters;
  const int rank = forest->mpirank, size = forest->mpisize;
  const t8_gloidx_t *offsets = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const t8_gloidx_t *tree_offsets = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
  const t8_linearidx_t *first_descs = (const t8_linearidx_t *) t8_shmem_array_get_array (forest->global_first_desc);
  const auto has_elements = [&] (const int r) { return offsets[r] < offsets[r + 1]; };
  /* The first rank in [from, to) with elements, or to */
  const auto first_with_elements = [&] (int from, const int to) {
    while (from < to && !has_elements (from)) {
      from++;
    }
    return from;
  };
  /* The first descendants of a process with deleted elements may lie in a later tree than its first tree */
  const bool check_positions = !forest->incomplete_trees;
  const auto position = [&] (const int r) {
    return t8_forest_position { t8_offset_first (r, tree_offsets), first_descs[r] };
  };

  bool match = splitters.rank == rank && splitters.first_element == offsets[rank]
               && splitters.global_num_elements == offsets[size] && splitters.has_elements == has_elements (rank);
  if (check_positions) {
    const int next = first_with_elements (rank + 1, size);
    match = match && (!splitters.has_elements || splitters.first_position == position (rank))
            && splitters.next_position == (next < size ? position (next) : splitters.end_position);
  }
  size_t ilevel = 0;
  for (int lo = 0, hi = size; hi - lo > 1; ilevel++) {
    const int mid = lo + (hi - lo + 1) / 2;
    if (ilevel >= splitters.levels.size ()) {
      return 0;
    }
    const t8_forest_bisection &level = splitters.levels[ilevel];
    const int first_upper = first_with_elements (mid, hi);
    match = match && level.lo == lo && level.mid == mid && level.hi == hi && level.mid_element == offsets[mid]
            && level.lower_has_elements == (first_with_elements (lo, mid) < mid)
            && level.upper_has_elements == (first_upper < hi);
    if (check_positions && first_upper < hi) {
      match = match && level.mid_position == position (first_upper);
    }
    (rank < mid ? hi : lo) = mid;
  }
  return match && ilevel == splitters.levels.size ();
}

void
t8_forest_element_positions (const t8_forest_t forest, const t8_gloidx_t gtreeid, const t8_eclass_t eclass,
                             const t8_element_t *element, t8_forest_position *first, t8_forest_position *last)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_element_t *desc;
  scheme->element_new (eclass, 1, &desc);
  scheme->element_get_first_descendant (eclass, element, desc, forest->maxlevel);
  *first = { gtreeid, scheme->element_get_linear_id (eclass, desc, forest->maxlevel) };
  scheme->element_get_last_descendant (eclass, element, desc, forest->maxlevel);
  *last = { gtreeid, scheme->element_get_linear_id (eclass, desc, forest->maxlevel) };
  scheme->element_destroy (eclass, 1, &desc);
}

int
t8_forest_route_halves_positions (const t8_forest_bisection &level, const t8_forest_position &first,
                                  const t8_forest_position &last, const bool unique)
{
  T8_ASSERT (first <= last);
  /* The lower half owns all positions of the interval before the first position of the upper half, or all of them
   * if the upper half has no elements. */
  int halves = T8_FOREST_ROUTE_NONE;
  if (level.lower_has_elements && (!level.upper_has_elements || first < level.mid_position)) {
    halves |= T8_FOREST_ROUTE_LOWER;
  }
  if (level.upper_has_elements && level.mid_position <= last) {
    halves |= T8_FOREST_ROUTE_UPPER;
  }
  return unique && halves == T8_FOREST_ROUTE_BOTH ? T8_FOREST_ROUTE_NONE : halves;
}

int
t8_forest_route_halves_elements (const t8_forest_bisection &level, const t8_gloidx_t first, const t8_gloidx_t last,
                                 const bool unique)
{
  T8_ASSERT (first <= last);
  int halves = T8_FOREST_ROUTE_NONE;
  if (level.lower_has_elements && first < level.mid_element) {
    halves |= T8_FOREST_ROUTE_LOWER;
  }
  if (level.upper_has_elements && level.mid_element <= last) {
    halves |= T8_FOREST_ROUTE_UPPER;
  }
  return unique && halves == T8_FOREST_ROUTE_BOTH ? T8_FOREST_ROUTE_NONE : halves;
}

void
t8_forest_route_exchange (const t8_forest_t forest, const t8_forest_bisection &level, const std::vector<char> &send,
                          std::vector<char> &recv)
{
  const int rank = forest->mpirank;
  const bool is_lower = rank < level.mid;
  const int j = rank - (is_lower ? level.lo : level.mid);
  const int num_upper = level.hi - level.mid;
  int target, sources[2], num_sources = 0;
  if (is_lower) {
    /* The last rank of a lower half with one rank more sends its items to the last rank of the upper half */
    target = j < num_upper ? level.mid + j : level.hi - 1;
    if (j < num_upper) {
      sources[num_sources++] = target;
    }
  }
  else {
    target = sources[num_sources++] = level.lo + j;
    if (j == num_upper - 1 && level.mid - level.lo > num_upper) {
      sources[num_sources++] = level.mid - 1;
    }
  }
  T8_ASSERT (send.size () <= (size_t) INT_MAX);
  sc_MPI_Request request;
  int mpiret = sc_MPI_Isend (const_cast<char *> (send.data ()), static_cast<int> (send.size ()), sc_MPI_BYTE, target,
                             T8_MPI_FOREST_ROUTE, forest->mpicomm, &request);
  SC_CHECK_MPI (mpiret);
  for (int isource = 0; isource < num_sources; isource++) {
    sc_MPI_Status status;
    mpiret = sc_MPI_Probe (sources[isource], T8_MPI_FOREST_ROUTE, forest->mpicomm, &status);
    SC_CHECK_MPI (mpiret);
    int num_bytes;
    mpiret = sc_MPI_Get_count (&status, sc_MPI_BYTE, &num_bytes);
    SC_CHECK_MPI (mpiret);
    const size_t offset = recv.size ();
    recv.resize (offset + num_bytes);
    mpiret = sc_MPI_Recv (recv.data () + offset, num_bytes, sc_MPI_BYTE, sources[isource], T8_MPI_FOREST_ROUTE,
                          forest->mpicomm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Waitall (1, &request, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
}
