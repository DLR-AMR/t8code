#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/levelmultiindex.hxx"
#include <t8.h>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace t8_mra::adapt
{

/// Reset all per-pass multiscale state.
template <typename MS>
void
clear_state (MS &mra)
{
  mra.d_map.erase_all ();
  mra.td_set.erase_all ();
  mra.refinement_set.erase_all ();
  mra.coarsening_set.erase_all ();
}

/// Number of leaves marked for refinement in [min_level, max_level).
template <typename MS>
unsigned int
num_refinement_marks (MS &mra, int min_level, int max_level)
{
  auto num = 0u;
  for (auto l = min_level; l < max_level; ++l)
    num += mra.refinement_set[l].size ();
  return num;
}

/// Reconstruct children data (inverse two-scale, zero details) for the marks in
/// refinement_set, then realize them with one forest adapt.
template <typename MS>
void
apply_refinement (MS &mra, int min_level, int max_level, int recursive)
{
  mra.d_map.erase_all ();
  for (auto l = min_level; l < max_level; ++l)
    for (const auto &lmi : mra.refinement_set[l])
      mra.d_map.insert (lmi, typename MS::detail_t {});

  mra.inverse_multiscale_transformation (min_level, max_level);
  mra.grid.adapt (MS::static_refinement_callback, recursive);
}

/// Stand-in for callers that realize their marks each round (balance) and never
/// descend a prior-refinements path.
struct no_prior_marks
{
  template <typename Lmi>
  bool
  contains (const Lmi &) const
  {
    return false;
  }
};

/**
 * @brief Find the covering leaf of a same-level neighbour and refine it
 *
 * Walks up to the covering leaf, then descends towards the neighbour through
 * prior_refinements (marks treated as performed). A covering leaf more than
 * max_level_gap levels coarser than the neighbour is refined one level (0:
 * grading, exact match; 1: 2:1 balance).
 *
 * @return 1 on a new mark, 0 if nothing to do, -1 if no covering leaf is local
 */
template <typename MS, typename Lmi, typename PriorRefinements>
int
refine_covering_leaf (MS &mra, const Lmi &neigh_lmi, int min_level, unsigned int max_level_gap,
                      const PriorRefinements &prior_refinements)
{
  auto *lmi_map = mra.get_lmi_map ();

  auto walk = neigh_lmi;
  while (walk.level () > 0 && !lmi_map->contains (walk))
    walk = t8_mra::parent_lmi (walk);

  if (!lmi_map->contains (walk))
    return -1;

  while (walk.level () + max_level_gap < neigh_lmi.level () && prior_refinements.contains (walk)) {
    auto down = neigh_lmi;
    while (down.level () > walk.level () + 1)
      down = t8_mra::parent_lmi (down);
    walk = down;
  }

  if (walk.level () + max_level_gap < neigh_lmi.level () && static_cast<int> (walk.level ()) >= min_level
      && !mra.refinement_set.contains (walk)) {
    mra.refinement_set.insert (walk);
    return 1;
  }

  return 0;
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
template <typename MS, typename PriorRefinements>
unsigned int
exchange_refine_requests (MS &mra, const std::vector<std::vector<size_t>> &outgoing, int min_level,
                          unsigned int max_level_gap, const PriorRefinements &prior_refinements)
{
  using levelmultiindex = typename MS::levelmultiindex;

  int mpisize;
  sc_MPI_Comm_size (mra.grid.comm, &mpisize);

  std::vector<int> send_counts (mpisize);
  std::transform (outgoing.begin (), outgoing.end (), send_counts.begin (),
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
    for (const auto index : batch)
      if (refine_covering_leaf (mra, levelmultiindex (index), min_level, max_level_gap, prior_refinements) > 0)
        ++num_new_marks;

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
template <typename MS, typename PriorRefinements>
unsigned int
neighbour_prediction (MS &mra, int min_level, const PriorRefinements &prior_refinements)
{
  int mpirank, mpisize;
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
