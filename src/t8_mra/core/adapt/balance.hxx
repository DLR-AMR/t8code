#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/adapt/grading.hxx"
#include <t8.h>

#include <vector>

namespace t8_mra::adapt
{

/**
 * @brief One balancing round: refine covering leaves across faces
 *
 * Every leaf resolves its face neighbours to their covering leaf; a covering
 * leaf more than one level coarser is refined one level (children data =
 * inverse two-scale with zero details, so the data is unchanged).
 *
 * @return Number of leaves marked in this round
 */
template <typename MS>
unsigned int
balance_round (MS &mra)
{
  clear_state (mra);

  int mpirank, mpisize;
  sc_MPI_Comm_rank (mra.grid.comm, &mpirank);
  sc_MPI_Comm_size (mra.grid.comm, &mpisize);
  std::vector<std::vector<size_t>> outgoing (mpisize);

  mra.grid.for_each_face_neigh (
    [] (const auto &lmi) { return lmi.level () >= 2; },
    [&] (const auto &, t8_eclass_t tree_class, t8_gloidx_t neigh_gtreeid, t8_element_t *neigh_element,
         const auto &neigh_lmi) {
      if (refine_covering_leaf (mra, neigh_lmi, 0, 1u, no_prior_marks {}) < 0 && mpisize > 1) {
        const auto owner = mra.grid.find_owner (neigh_gtreeid, neigh_element, tree_class);
        if (owner != mpirank)
          outgoing[owner].push_back (neigh_lmi.index);
      }
    });

  if (mpisize > 1)
    exchange_refine_requests (mra, outgoing, 0, 1u, no_prior_marks {});

  const auto num_marked = mra.grid.global_num_marks (num_refinement_marks (mra, 0, mra.grid.maximum_level));
  if (num_marked == 0) {
    clear_state (mra);
    return 0;
  }

  apply_refinement (mra, 0, mra.grid.maximum_level, 0);
  clear_state (mra);

  return num_marked;
}

/**
 * @brief Restore the 2:1 face balance of the grid
 *
 * Rounds iterate until no leaf has a face neighbour more than one level coarser;
 * a jump of k levels resolves in k-1 rounds. Terminates: every round refines at
 * least one leaf and levels are bounded by max_level.
 */
template <typename MS>
void
balance (MS &mra)
{
  auto rounds = 0;
  while (balance_round (mra) > 0)
    t8_debugf ("MRA balance round %d\n", rounds++);

  if (rounds > 0)
    mra.grid.repartition ();
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
