#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>

#include <t8.h>

#include "t8_mra/core/adapt/grading.hxx"
#include "t8_mra/criteria/coarsening_criterion.hxx"
#include "t8_mra/criteria/refinement_criterion.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

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
template <typename TMultiscale, typename TCriterion = harten_prediction>
  requires refinement_criterion<TCriterion, TMultiscale>
void
refine (TMultiscale &mra, int min_level, int max_level, TCriterion criterion = {})
{
  if constexpr (criterion_has_prepare<TCriterion, TMultiscale>)
    criterion.prepare (mra);

  clear_state (mra);

  mra.multiscale_transformation (0, max_level);

  auto num_families = 0u;
  for (auto l = 0; l < max_level; ++l) {
    for (const auto &[lmi, _] : mra.d_map[l]) {
      ++num_families;

      const auto flags = criterion (mra, lmi);

      if (flags.grade_neighbours)
        mra.td_set.insert (lmi);

      if (l + 1 >= min_level && l < max_level - 1 && flags.refine_children)
        for (const auto &child : t8_mra::children_lmi (lmi))
          mra.refinement_set.insert (child);
    }
  }

  typename TMultiscale::level_set prior_refinements (mra.maximum_level ());

  for (auto round = 0;; ++round) {
    const auto new_marks = neighbour_prediction (mra, min_level, prior_refinements);

    t8_debugf ("MRA refine grading round %d: %u new marks\n", round, new_marks);
    if (new_marks == 0)
      break;

    prior_refinements = mra.refinement_set;

    typename TMultiscale::index_set stopped;
    for (auto l = 0; l < max_level; ++l)
      for (const auto &lmi : mra.td_set[l]) {
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
