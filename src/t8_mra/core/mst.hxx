#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "t8_mra/core/shape/mst_policy.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/num/mask_coefficients.hxx"

namespace t8_mra
{

/**
 * @brief Two-scale (multiscale) transform operations on a levelindex_map.
 *
 *   - two_scale_family: per-family kernel (children -> parent + details).
 *   - multiscale_transformation: forward, non-destructive; details to d_map.
 *   - multiscale_decomposition: forward, destructive; collapse to l_min.
 *   - inverse_multiscale_transformation: reconstruct children from details.
 *
 * Element-specific behaviour is routed through the ordering and scaling policies.
 */
template <typename TElement, typename TDetail = detail_data<TElement::Shape, TElement::U_DIM, TElement::P_DIM>,
          typename TOrderingPolicy = ordering_policy<TElement::Shape>,
          typename TScalingPolicy = mst_scaling_policy<TElement::Shape>>
class mst {
 public:
  using element_t = TElement;
  using detail_t = TDetail;
  using levelmultiindex = t8_mra::levelmultiindex<TElement::Shape>;

  static constexpr auto Shape = TElement::Shape;
  static constexpr unsigned int U_DIM = TElement::U_DIM;
  static constexpr unsigned int DOF = TElement::DOF;

  using mask_t = two_scale_mask<Shape, TElement::P_DIM>;

  /// Two-scale mask coefficients, computed once from the reference basis.
  mask_t mask;

  mst ()
  {
    t8_mra::compute_mask<Shape, TElement::P_DIM> (mask);
  }

  void
  two_scale_family (const std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                    detail_t &data_on_coarse) const
  {
    two_scale_family (data_on_siblings, data_on_coarse, mask);
  }

  void
  inverse_two_scale_family (const element_t &data_on_coarse, const detail_t &details,
                            std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings) const
  {
    inverse_two_scale_family (data_on_coarse, details, data_on_siblings, mask);
  }

  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max,
                             levelindex_map<levelmultiindex, element_t> &lmi_map,
                             levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    multiscale_transformation (l_min, l_max, lmi_map, d_map, mask);
  }

  void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> &lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    multiscale_decomposition (l_min, l_max, lmi_map, d_map, mask);
  }

  template <typename TKeep, typename TCollapsed>
  void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> &lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map, TKeep &&keep, TCollapsed &&collapsed) const
  {
    multiscale_decomposition (l_min, l_max, lmi_map, d_map, mask, std::forward<TKeep> (keep),
                              std::forward<TCollapsed> (collapsed));
  }

  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max,
                                     levelindex_map<levelmultiindex, element_t> &lmi_map,
                                     levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    inverse_multiscale_transformation (l_min, l_max, lmi_map, d_map, mask);
  }

  /**
   * @brief Parents of the families present on a level, one entry each.
   *
   * Only the leaf that is child 0 of its family contributes, so no deduplication is
   * needed; a family whose child 0 is not a leaf here is incomplete and would be
   * skipped anyway. The result is a snapshot, which the destructive sweep needs.
   */
  static void
  collect_parents (const levelindex_map<levelmultiindex, element_t> &lmi_map, unsigned int level,
                   std::vector<levelmultiindex> &parents)
  {
    parents.clear ();
    parents.reserve (lmi_map.size (level) / levelmultiindex::NUM_CHILDREN);

    for (const auto &[lmi, _] : lmi_map[level])
      if (lmi.child_id (level) == 0)
        parents.push_back (t8_mra::parent_lmi (lmi));
  }

  /**
   * @brief Two-scale transform of one complete family (children -> parent + details).
   *
   *   u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
   *   d[k][i]     = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
   *
   * @param data_on_siblings The NUM_CHILDREN children of the family.
   * @param data_on_coarse   Output parent: u_coeffs, d_coeffs, vol, order.
   * @param mask_coefficients Two-scale mask matrices M[k].
   */
  static void
  two_scale_family (const std::array<const element_t *, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                    detail_t &data_on_coarse, const mask_t &mask_coefficients)
  {
    const double scaling_factor = TScalingPolicy::forward_scaling_factor (levelmultiindex::NUM_CHILDREN);

    for (auto u = 0u; u < U_DIM; ++u) {
      std::array<double, DOF> u_parent;

      // Parent coefficients: u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
      for (auto i = 0u; i < DOF; ++i) {
        auto sum = 0.0;

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          const auto &Mk_column = mask_coefficients.transposed[k][i];
          const auto &uk = data_on_siblings[k]->u_coeffs;

          for (auto j = 0u; j < DOF; ++j)
            sum += uk[element_t::dg_idx (u, j)] * Mk_column[j];
        }

        u_parent[i] = sum * scaling_factor;
        data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = u_parent[i];
      }

      // Detail coefficients: d[k][i] = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        const auto &Mk = mask_coefficients.m[k];
        const auto &uk = data_on_siblings[k]->u_coeffs;

        for (auto i = 0u; i < DOF; ++i) {
          const auto &Mk_row = Mk[i];
          auto sum = 0.0;

          for (auto j = 0u; j < DOF; ++j)
            sum += Mk_row[j] * u_parent[j];

          data_on_coarse.d_coeffs[detail_t::wavelet_idx (k, u, i)] = uk[element_t::dg_idx (u, i)] - sum;
        }
      }
    }

    data_on_coarse.vol = data_on_siblings[0]->vol * levelmultiindex::NUM_CHILDREN;
    data_on_coarse.order = data_on_siblings[0]->order;

    TOrderingPolicy::adjust_parent_order (data_on_coarse);
  }

  /// Convenience form for callers that hold the siblings by value.
  static void
  two_scale_family (const std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                    detail_t &data_on_coarse, const mask_t &mask_coefficients)
  {
    std::array<const element_t *, levelmultiindex::NUM_CHILDREN> siblings;

    for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
      siblings[k] = &data_on_siblings[k];

    two_scale_family (siblings, data_on_coarse, mask_coefficients);
  }

  /**
   * @brief Inverse two-scale transform of one family (parent + details -> children).
   *
   *   u_child[k][i] = d[k][i] + Σ_j M[k](i,j) * u_parent[j]
   *
   * @param data_on_coarse   The parent cell.
   * @param details          The family's detail coefficients.
   * @param data_on_siblings Output children.
   * @param mask_coefficients Two-scale mask matrices M[k].
   */
  static void
  inverse_two_scale_family (const element_t &data_on_coarse, const detail_t &details,
                            std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                            const mask_t &mask_coefficients)
  {
    const double inv_scaling_factor = TScalingPolicy::inverse_scaling_factor ();
    const auto &u_parent = data_on_coarse.u_coeffs;

    for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
      const auto &Mk = mask_coefficients.m[k];
      auto &child = data_on_siblings[k];

      for (auto u = 0u; u < U_DIM; ++u) {
        for (auto i = 0u; i < DOF; ++i) {
          const auto &Mk_row = Mk[i];
          auto sum = 0.0;

          for (auto j = 0u; j < DOF; ++j)
            sum += u_parent[element_t::dg_idx (u, j)] * Mk_row[j];

          child.u_coeffs[element_t::dg_idx (u, i)]
            = details.d_coeffs[detail_t::wavelet_idx (k, u, i)] + sum * inv_scaling_factor;
        }
      }

      child.vol = data_on_coarse.vol / levelmultiindex::NUM_CHILDREN;
      TOrderingPolicy::adjust_child_order (child, k, data_on_coarse);
    }
  }

  /**
   * @brief Volume-scaled detail 2-norm per component.
   *
   * @param  detail A family's detail coefficients and volume.
   * @return Per-component detail norm.
   */
  [[nodiscard]] static std::array<double, U_DIM>
  detail_norm (const detail_t &detail)
  {
    std::array<double, U_DIM> norm = {};
    const auto &details = detail.d_coeffs;

    for (auto u = 0u; u < U_DIM; ++u) {
      auto norm_sq = 0.0;

      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
        for (auto i = 0u; i < DOF; ++i) {
          const auto d = details[detail_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }

      norm[u] = std::sqrt (norm_sq * TScalingPolicy::detail_norm_scale (detail.vol));
    }

    return norm;
  }

  /**
   * @brief Forward transform of every complete family in (l_min, l_max]; writes
   *        their details to d_map, leaves lmi_map unchanged.
   *
   * @param l_min, l_max     Level range, exclusive of l_min.
   * @param lmi_map          Single-scale leaves.
   * @param d_map            Output details per family.
   * @param mask_coefficients Two-scale mask matrices M[k].
   */
  static void
  multiscale_transformation (unsigned int l_min, unsigned int l_max,
                             levelindex_map<levelmultiindex, element_t> &lmi_map,
                             levelindex_map<levelmultiindex, detail_t> &d_map,
                             const mask_t &mask_coefficients)
  {
    std::vector<levelmultiindex> parents;
    detail_t data_on_coarse;
    std::array<const element_t *, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_max; l > l_min; --l) {
      d_map[l - 1].reserve (lmi_map.size (l) / levelmultiindex::NUM_CHILDREN);

      collect_parents (lmi_map, l, parents);

      for (const auto &lmi : parents) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // Incomplete families (siblings on finer levels) carry no detail
        // information.
        auto family_complete = true;
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          data_on_siblings[k] = lmi_map.find (siblings_lmi[k]);

          if (data_on_siblings[k] == nullptr) {
            family_complete = false;
            break;
          }
        }

        if (!family_complete)
          continue;

        two_scale_family (data_on_siblings, data_on_coarse, mask_coefficients);
        d_map.insert (lmi, data_on_coarse);
      }
    }
  }

  /**
   * @brief Collapse each complete family down to l_min: replace it by its parent
   *        in lmi_map (erase children) and write its details to d_map.
   *
   * @param l_min, l_max     Level range, exclusive of l_min.
   * @param lmi_map          Leaves; collapsed in place.
   * @param d_map            Output details per family.
   * @param mask_coefficients Two-scale mask matrices M[k].
   * @param keep             keep(parent lmi) after its detail is in d_map: true
   *                         leaves the family refined (thresholded coarsening).
   * @param collapsed        collapsed(child lmi) per consumed child.
   */
  template <typename TKeep, typename TCollapsed>
  static void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> &lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map,
                            const mask_t &mask_coefficients, TKeep &&keep, TCollapsed &&collapsed)
  {
    std::vector<levelmultiindex> parents;
    detail_t data_on_coarse;
    std::array<const element_t *, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_max; l > l_min; --l) {
      d_map[l - 1].reserve (lmi_map.size (l) / levelmultiindex::NUM_CHILDREN);

      collect_parents (lmi_map, l, parents);

      for (const auto &lmi : parents) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // On an adaptive grid a family may be incomplete: some siblings stayed
        // refined on finer levels. Such families cannot be two-scale transformed.
        auto family_complete = true;
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          data_on_siblings[k] = lmi_map.find (siblings_lmi[k]);

          if (data_on_siblings[k] == nullptr) {
            family_complete = false;
            break;
          }
        }

        if (!family_complete)
          continue;

        two_scale_family (data_on_siblings, data_on_coarse, mask_coefficients);
        d_map.insert (lmi, data_on_coarse);

        if (keep (lmi))
          continue;

        // The lmi_map leaf keeps only single-scale data (slice off d_coeffs).
        lmi_map.insert (lmi, static_cast<const element_t &> (data_on_coarse));

        // Consume only this family's children; members of skipped (incomplete)
        // families must stay in the map as leaves.
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          lmi_map.erase (siblings_lmi[k]);
          collapsed (siblings_lmi[k]);
        }
      }
    }
  }

  /** @brief Unconditional decomposition: every complete family collapses. */
  static void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> &lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map,
                            const mask_t &mask_coefficients)
  {
    multiscale_decomposition (
      l_min, l_max, lmi_map, d_map, mask_coefficients, [] (const auto & /*unused*/) { return false; },
      [] (const auto & /*unused*/) {});
  }

  /**
   * @brief Reconstruct children from parent and detail coefficients over
   *        [l_min, l_max); moves the data from d_map back into lmi_map.
   *
   *   u_child[k][i] = d[k][i] + Σ_j M[k](i,j) * u_parent[j]
   *
   * @param l_min, l_max     Level range, exclusive of l_max.
   * @param lmi_map          Parents in, children out.
   * @param d_map            Details in; consumed.
   * @param mask_coefficients Two-scale mask matrices M[k].
   */
  static void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max,
                                     levelindex_map<levelmultiindex, element_t> &lmi_map,
                                     levelindex_map<levelmultiindex, detail_t> &d_map,
                                     const mask_t &mask_coefficients)
  {
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_min; l < l_max; ++l) {
      lmi_map[l + 1].reserve (d_map[l].size ());

      for (const auto &[lmi, d] : d_map[l]) {
        const auto children_lmi = t8_mra::children_lmi (lmi);

        inverse_two_scale_family (lmi_map.get (lmi), d, data_on_siblings, mask_coefficients);

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          lmi_map.insert (children_lmi[k], data_on_siblings[k]);

        lmi_map.erase (lmi);
      }

      d_map.erase (l);
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
