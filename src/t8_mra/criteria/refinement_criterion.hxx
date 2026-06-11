#pragma once

#ifdef T8_ENABLE_MRA

#include <cmath>
#include <concepts>

namespace t8_mra
{

/**
 * @brief Requirements for a refinement criterion
 *
 * A refinement criterion decides per leaf family (identified by the parent
 * lmi, whose detail data is available in mra.d_map) how the grid refines:
 *
 *   - refine_neighbours(mra, lmi): the grid around the family is graded so
 *     that all face neighbours reach the family's leaf level.
 *   - refine(mra, lmi): the family's children are refined one further
 *     level.
 *
 * Optionally a criterion can provide prepare(mra), which is called once at
 * the beginning of every refine() call.
 */
template <typename C, typename MRA>
concept refinement_criterion = requires (C c, MRA &mra, const typename MRA::levelmultiindex &lmi) {
  { c.refine_neighbours (mra, lmi) } -> std::convertible_to<bool>;
  { c.refine (mra, lmi) } -> std::convertible_to<bool>;
};

/**
 * @brief Example refinement criterion: Harten's prediction
 *
 * Significant details (hard threshold) grade their neighbourhood, a
 * steep-gradient detail additionally refines the family's children:
 *
 *   refine_neighbours:  max_u ||d_u|| / c_scaling_u  >           c_thresh * eps(lmi)
 *   refine:             max_u ||d_u|| / c_scaling_u  >  2^(P+1) * c_thresh * eps(lmi)
 *
 * prepare() computes the global scaling factors c_scaling (eq. 2.39).
 */
struct harten_prediction
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
  refine_neighbours (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    return mra.scaled_detail_norm (lmi) > c_thresh * mra.local_threshold_value (lmi, gamma);
  }

  template <typename MRA>
  bool
  refine (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    const auto steep_factor = std::pow (2.0, static_cast<int> (MRA::P_DIM) + 1);
    return mra.scaled_detail_norm (lmi) > steep_factor * c_thresh * mra.local_threshold_value (lmi, gamma);
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
