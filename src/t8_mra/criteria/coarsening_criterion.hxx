#pragma once

#ifdef T8_ENABLE_MRA

#include <concepts>

namespace t8_mra
{

/**
 * @brief Requirements for a coarsening criterion
 *
 * A coarsening criterion decides per leaf family (identified by the parent
 * lmi, whose detail data is available in mra.d_map) whether the family's
 * detail information is essential:
 *
 *   - significant(mra, lmi) == true:  the family keeps its details and
 *     stays refined.
 *   - significant(mra, lmi) == false: the details are discarded and the
 *     family is coarsened into its parent.
 *
 * Optionally a criterion can provide prepare(mra), which is called once at
 * the beginning of every coarsen() call (e.g. to compute global
 * normalization factors).
 */
template <typename TCriterion, typename TMultiscale>
concept coarsening_criterion
  = requires (TCriterion criterion, TMultiscale &mra, const typename TMultiscale::levelmultiindex &lmi) {
      { criterion.significant (mra, lmi) } -> std::convertible_to<bool>;
    };

/**
 * @brief Detect optional prepare() hook of a criterion
 */
template <typename TCriterion, typename TMultiscale>
concept criterion_has_prepare = requires (TCriterion criterion, TMultiscale &mra) { criterion.prepare (mra); };

/**
 * @brief Default coarsening criterion: hard thresholding
 *
 * Uses the level-dependent threshold of Veli eq. (2.44) on the scaled
 * detail norms:
 *
 *   significant:  max_u ||d_u|| / c_scaling_u  >  c_thresh * eps(lmi)
 *
 * prepare() computes the global scaling factors c_scaling (eq. 2.39).
 */
struct hard_thresholding
{
  /// Threshold constant
  double c_thresh = 1.0;
  /// Expected order of convergence
  int gamma = 1;

  template <typename TMultiscale>
  void
  prepare (TMultiscale &mra)
  {
    mra.c_scaling = mra.threshold_scaling_factor ();
  }

  template <typename TMultiscale>
  [[nodiscard]] bool
  significant (TMultiscale &mra, const typename TMultiscale::levelmultiindex &lmi)
  {
    return mra.scaled_detail_norm (lmi) > c_thresh * mra.local_threshold_value (lmi, gamma);
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
