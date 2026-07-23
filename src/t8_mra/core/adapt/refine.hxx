#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/adapt/grading.hxx"
#include "t8_mra/criteria/refinement_criterion.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include <t8.h>

#include <algorithm>

namespace t8_mra::adapt
{

/**
 * @brief Adaptive refinement from min_level up to max_level
 *
 * The non-destructive transform yields every family's details; the criterion
 * grades neighbourhoods (td_set) and refines steep families' children. One pass
 * suffices (new children carry zero details); a grading fixpoint then pulls
 * covering leaves up one level per round. All marks are realized by one recursive
 * forest adapt (children data = inverse two-scale with zero details).
 */
template <typename MS, typename Criterion = harten_prediction>
  requires refinement_criterion<Criterion, MS>
void
refine (MS &mra, int min_level, int max_level, Criterion criterion = {})
{
  if constexpr (criterion_has_prepare<Criterion, MS>)
    criterion.prepare (mra);

  clear_state (mra);

  mra.multiscale_transformation (0, max_level);

  auto num_families = 0u;
  for (auto L = 0; L < max_level; ++L) {
    for (const auto &[lmi, _] : mra.d_map[L]) {
      ++num_families;

      const auto flags = criterion (mra, lmi);

      if (flags.grade_neighbours)
        mra.td_set.insert (lmi);

      if (L < max_level - 1 && flags.refine_children)
        for (const auto &child : t8_mra::children_lmi (lmi))
          mra.refinement_set.insert (child);
    }
  }

  for (auto l = 0; l < min_level; ++l)
    mra.refinement_set.erase (l);

  auto prior_refinements = mra.refinement_set;
  prior_refinements.erase_all ();

  for (auto round = 0;; ++round) {
    const auto new_marks = neighbour_prediction (mra, min_level, prior_refinements);
    t8_debugf ("MRA refine grading round %d: %u new marks\n", round, new_marks);
    if (new_marks == 0)
      break;

    prior_refinements = mra.refinement_set;

    typename MS::index_set stopped;
    for (auto L = 0; L < max_level; ++L)
      for (const auto &lmi : mra.td_set[L]) {
        const auto children = t8_mra::children_lmi (lmi);
        if (std::any_of (children.begin (), children.end (),
                         [&] (const auto &child) { return mra.refinement_set.contains (child); }))
          stopped.insert (lmi);
      }
    for (const auto &lmi : stopped)
      mra.td_set.erase (lmi);
  }

  const auto num_marked = num_refinement_marks (mra, min_level, max_level);
  t8_debugf ("MRA refine analysis: %u leaf families, %u leaves marked\n", num_families, num_marked);

  if (mra.grid.global_num_marks (num_marked) == 0) {
    clear_state (mra);
    return;
  }

  apply_refinement (mra, min_level, max_level, 1);
  mra.grid.repartition ();

  clear_state (mra);
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
