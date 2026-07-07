#pragma once

#ifdef T8_ENABLE_MRA

#include <cmath>
#include <concepts>

namespace t8_mra
{

/**
 * @brief Per-family refinement decision
 *
 * The two flags are independent:
 *   - grade_neighbours: grade the surrounding grid so the family's face
 *     neighbours reach its leaf level.
 *   - refine_children: refine the family's children one further level.
 */
struct refinement_flags
{
  bool grade_neighbours;
  bool refine_children;
};

/**
 * @brief Requirements for a refinement criterion
 *
 * operator()(mra, lmi) returns the refinement_flags for a leaf family
 * (identified by the parent lmi, whose detail is in mra.d_map). Optionally a
 * criterion provides prepare(mra), called once at the start of every refine().
 */
template <typename C, typename MRA>
concept refinement_criterion = requires (C c, MRA &mra, const typename MRA::levelmultiindex &lmi) {
  { c (mra, lmi) } -> std::convertible_to<refinement_flags>;
};

/**
 * @brief Example refinement criterion: Harten's prediction
 *
 * On the scaled detail norm N = max_u ||d_u|| / c_scaling_u and the level
 * threshold eps(lmi):
 *   grade_neighbours:  N >            c_thresh * eps
 *   refine_children:   N >  2^(P+1) * c_thresh * eps
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
  refinement_flags
  operator() (MRA &mra, const typename MRA::levelmultiindex &lmi)
  {
    const auto norm = mra.scaled_detail_norm (lmi);
    const auto threshold = c_thresh * mra.local_threshold_value (lmi, gamma);
    const auto steep_factor = std::pow (2.0, static_cast<int> (MRA::P_DIM) + 1);

    return { norm > threshold, norm > steep_factor * threshold };
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
