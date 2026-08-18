#pragma once

#ifdef T8_ENABLE_MRA

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
 * grades neighbourhoods (td_set) and refines steep families' children. One
 * grading pass marks the full refinement path of every coarser covering
 * leaf around td_set; all marks are realized by one recursive forest adapt
 * (children data = inverse two-scale with zero details).
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

  grade_neighbours (mra, min_level, 0u, [&] (const auto &lmi) {
    return lmi.level () != 0 && mra.td_set.contains (t8_mra::parent_lmi (lmi));
  });

  const auto num_marked = static_cast<unsigned int> (mra.refinement_set.size ());
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
