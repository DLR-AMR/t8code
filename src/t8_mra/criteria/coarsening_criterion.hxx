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
template <typename C, typename MRA>
concept coarsening_criterion = requires (C c, MRA &mra, const typename MRA::levelmultiindex &lmi) {
  { c.significant (mra, lmi) } -> std::convertible_to<bool>;
};

/**
 * @brief Detect optional prepare() hook of a criterion
 */
template <typename C, typename MRA>
concept criterion_has_prepare = requires (C c, MRA &mra) { c.prepare (mra); };

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
  /// Expected order of convergence (enters the level-dependent threshold)
  int gamma = 1;

  template <typename MRA>
  void
  prepare (MRA &mra)
  {
    /// Scaling due to (2.39)
    mra.c_scaling = mra.threshold_scaling_factor ();
  }

  template <typename MRA>
  bool
  significant (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    return mra.scaled_detail_norm (lmi) > c_thresh * mra.local_threshold_value (lmi, gamma);
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
