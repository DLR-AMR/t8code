#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <vector>

#include "sc_mpi.h"

#include "t8.h"
#include "t8_eclass/t8_eclass.h"
#include "t8_element/t8_element.h"

#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra::adapt
{

/// Reset all per-pass multiscale state.
template <typename TMultiscale>
void
clear_state (TMultiscale &mra)
{
  mra.d_map.erase_all ();
  mra.td_set.erase_all ();
  mra.refinement_set.erase_all ();
  mra.coarsening_set.erase_all ();
}

/// Number of leaves marked for refinement in [min_level, max_level).
template <typename TMultiscale>
[[nodiscard]] unsigned int
num_refinement_marks (TMultiscale &mra, int min_level, int max_level)
{
  auto num = 0u;

  for (auto l = min_level; l < max_level; ++l)
    num += mra.refinement_set[l].size ();

  return num;
}

/// Reconstruct children data (inverse two-scale, zero details) for the marks in
/// refinement_set, then realize them with one forest adapt.
template <typename TMultiscale>
void
apply_refinement (TMultiscale &mra, int min_level, int max_level, int recursive)
{
  mra.d_map.erase_all ();
  for (auto l = min_level; l < max_level; ++l)
    for (const auto &lmi : mra.refinement_set[l])
      mra.d_map.insert (lmi, typename TMultiscale::detail_t {});

  mra.inverse_multiscale_transformation (min_level, max_level);
  mra.grid.adapt (TMultiscale::static_refinement_callback, recursive);
}

/// Stand-in for callers that realize their marks each round (balance) and never
/// descend a prior-refinements path.
struct no_prior_marks
{
  template <typename TLmi>
  bool
  contains (const TLmi & /*unused*/) const
  {
    return false;
  }
};

/**
 * @brief Mark the refinement path from a neighbour's covering leaf to its level
 *
 * Walks up to the covering leaf, then marks every cell on the path back down,
 * so the resulting leaf ends at most max_level_gap levels coarser than the
 * neighbour (0: grading, exact match; 1: 2:1 balance). The path realizes in
 * one recursive adapt.
 *
 * @return Number of new marks, or -1 if no covering leaf is local
 */
template <typename TMultiscale, typename TLmi>
[[nodiscard]] int
mark_refinement_path (TMultiscale &mra, const TLmi &neigh_lmi, int min_level, unsigned int max_level_gap)
{
  auto *lmi_map = mra.get_lmi_map ();

  std::array<TLmi, 1u << TLmi::LEVEL_BITS> ancestor;
  auto walk = neigh_lmi;
  ancestor[walk.level ()] = walk;

  while (walk.level () > 0 && !lmi_map->contains (walk)) {
    walk = t8_mra::parent_lmi (walk);
    ancestor[walk.level ()] = walk;
  }

  if (!lmi_map->contains (walk))
    return -1;

  if (static_cast<int> (walk.level ()) < min_level)
    return 0;

  auto num_marks = 0;
  for (auto level = walk.level (); level + max_level_gap < neigh_lmi.level (); ++level)
    if (!mra.refinement_set.contains (ancestor[level])) {
      mra.refinement_set.insert (ancestor[level]);
      ++num_marks;
    }

  return num_marks;
}

/**
 * @brief Ship pull-up requests to their owner ranks, resolve received ones
 *
 * A same-level neighbour whose covering leaf is not local is sent to the rank
 * owning that region; the owner resolves it against its own lmi_map and marks
 * its own refinement_set. Collective.
 *
 * @return Number of new LOCAL marks created by received requests
 */
template <typename TMultiscale, typename TPriorRefinements>
unsigned int
exchange_refine_requests (TMultiscale &mra, const std::vector<std::vector<size_t>> &outgoing, int min_level,
                          unsigned int max_level_gap, const TPriorRefinements &prior_refinements)
{
  using levelmultiindex = typename TMultiscale::levelmultiindex;

  int mpisize;
  sc_MPI_Comm_size (mra.grid.comm, &mpisize);

  std::vector<int> send_counts;
  send_counts.reserve (mpisize);

  std::ranges::transform (outgoing, std::back_inserter (send_counts),
                          [] (const auto &list) { return static_cast<int> (list.size ()); });

  std::vector<int> recv_counts (mpisize, 0);
  sc_MPI_Alltoall (send_counts.data (), 1, sc_MPI_INT, recv_counts.data (), 1, sc_MPI_INT, mra.grid.comm);

  std::vector<std::vector<size_t>> incoming (mpisize);
  std::vector<sc_MPI_Request> requests;
  requests.reserve (2 * mpisize);

  for (auto rank = 0; rank < mpisize; ++rank) {
    if (recv_counts[rank] > 0) {
      incoming[rank].resize (recv_counts[rank]);
      requests.emplace_back ();

      sc_MPI_Irecv (incoming[rank].data (), recv_counts[rank] * sizeof (size_t), sc_MPI_BYTE, rank, 0, mra.grid.comm,
                    &requests.back ());
    }

    if (send_counts[rank] > 0) {
      requests.emplace_back ();
      sc_MPI_Isend (const_cast<size_t *> (outgoing[rank].data ()), send_counts[rank] * sizeof (size_t), sc_MPI_BYTE,
                    rank, 0, mra.grid.comm, &requests.back ());
    }
  }
  sc_MPI_Waitall (static_cast<int> (requests.size ()), requests.data (), sc_MPI_STATUSES_IGNORE);

  auto num_new_marks = 0u;
  for (const auto &batch : incoming)
    num_new_marks += static_cast<unsigned int> (std::ranges::count_if (batch, [&] (const auto &idx) {
      return refine_covering_leaf (mra, levelmultiindex (idx), min_level, max_level_gap, prior_refinements) > 0;
    }));

  return num_new_marks;
}

/**
 * @brief One grading round: same-level neighbours of marked families
 *
 * Every leaf of a family in td_set refines its coarser covering leaves by one
 * level; larger jumps resolve over repeated rounds against prior_refinements. A
 * neighbour whose covering leaf is remote is shipped to its owner. Collective:
 * the returned count is global.
 *
 * @return Global number of new marks in this round
 */
template <typename TMultiscale, typename TPriorRefinements>
[[nodiscard]] unsigned int
neighbour_prediction (TMultiscale &mra, int min_level, const TPriorRefinements &prior_refinements)
{
  int mpirank;
  int mpisize;
  sc_MPI_Comm_rank (mra.grid.comm, &mpirank);
  sc_MPI_Comm_size (mra.grid.comm, &mpisize);

  std::vector<std::vector<size_t>> outgoing (mpisize);
  auto num_new_marks = 0u;

  mra.grid.for_each_face_neigh (
    [&] (const auto &lmi) { return lmi.level () != 0 && mra.td_set.contains (t8_mra::parent_lmi (lmi)); },
    [&] (const auto &, t8_eclass_t tree_class, t8_gloidx_t neigh_gtreeid, t8_element_t *neigh_element,
         const auto &neigh_lmi) {
      const auto res = refine_covering_leaf (mra, neigh_lmi, min_level, 0u, prior_refinements);
      if (res > 0)
        ++num_new_marks;
      else if (res < 0 && mpisize > 1) {
        const auto owner = mra.grid.find_owner (neigh_gtreeid, neigh_element, tree_class);

        if (owner != mpirank)
          outgoing[owner].push_back (neigh_lmi.index);
      }
    });

  if (mpisize > 1)
    num_new_marks += exchange_refine_requests (mra, outgoing, min_level, 0u, prior_refinements);

  return mra.grid.global_num_marks (num_new_marks);
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
