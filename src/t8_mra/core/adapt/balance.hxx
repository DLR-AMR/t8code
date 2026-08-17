#pragma once

#ifdef T8_ENABLE_MRA

#include "t8.h"

#include "t8_mra/core/adapt/grading.hxx"

namespace t8_mra::adapt
{

/**
 * @brief Restore the 2:1 face balance of the grid
 *
 * Each round marks the refinement paths of all covering leaves more than one
 * level coarser than a face neighbour and realizes them (children data =
 * inverse two-scale with zero details, so the data is unchanged); rounds
 * iterate until the new leaves open no further violations. Terminates: every
 * round refines at least one leaf and levels are bounded by max_level.
 */
template <typename TMultiscale>
void
balance (TMultiscale &mra)
{
  auto rounds = 0;

  for (;; ++rounds) {
    clear_state (mra);

    const auto num_marked
      = mra.grid.global_num_marks (grade_neighbours (mra, 0, 1u, [] (const auto &lmi) { return lmi.level () >= 2; }));

    if (num_marked == 0)
      break;

    t8_debugf ("MRA balance round %d: %u marks\n", rounds, num_marked);
    apply_refinement (mra, 0, mra.grid.maximum_level, 1);
  }

  clear_state (mra);

  if (rounds > 0)
    mra.grid.repartition ();
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
