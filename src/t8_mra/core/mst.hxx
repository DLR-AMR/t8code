#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/shape/cartesian.hxx"
#include "t8_mra/core/shape/triangle.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelindex_set.hxx"
#include "t8_mra/num/mat.hxx"
#include "t8_mra/num/dg_basis.hxx"
#include "t8_mra/num/mask_coefficients.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace t8_mra
{

/**
 * @brief Multiscale transformation (MST) operations.
 *
 * The two-scale relation between a family of children and its parent is the
 * core of the method. Three operations build on it:
 *   - multiscale_transformation: forward analysis (parent + details into
 *     d_map), non-destructive. coarsen/refine run this to get the details
 *     their criteria act on.
 *   - inverse_multiscale_transformation: reconstruction (coarse + details ->
 *     children). refinement's realization step.
 *   - multiscale_decomposition: full destructive forward collapse to the
 *     coarsest level (single-scale -> multiscale representation).
 * two_scale_family is the per-family kernel they all share.
 *
 * Works for triangular and cartesian elements; element-specific behaviour is
 * routed through the ordering and scaling policies.
 */
template <typename TElement, typename TDetail = detail_data<TElement::Shape, TElement::U_DIM, TElement::P_DIM>,
          typename ordering_policy_t = ordering_policy<TElement::Shape>,
          typename scaling_policy_t = mst_scaling_policy<TElement::Shape>>
class mst {
 public:
  using element_t = TElement;  // leaf data (lmi_map)
  using detail_t = TDetail;    // leaf data + details (d_map)
  using levelmultiindex = t8_mra::levelmultiindex<TElement::Shape>;
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;

  static constexpr auto Shape = TElement::Shape;
  static constexpr unsigned int U_DIM = TElement::U_DIM;
  static constexpr unsigned int DOF = TElement::DOF;

  /// Two-scale mask coefficients, computed once from the reference basis.
  std::vector<t8_mra::mat> mask;

  mst ()
  {
    t8_mra::compute_mask<Shape, TElement::P_DIM> (mask);
  }

  /// Instance forwards baking the owned mask; the static overloads take an
  /// explicit mask (used by the standalone MST tests).
  void
  two_scale_family (const std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                    detail_t &data_on_coarse) const
  {
    two_scale_family (data_on_siblings, data_on_coarse, mask);
  }

  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max,
                             levelindex_map<levelmultiindex, element_t> *lmi_map,
                             levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    multiscale_transformation (l_min, l_max, lmi_map, d_map, mask);
  }

  void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> *lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    multiscale_decomposition (l_min, l_max, lmi_map, d_map, mask);
  }

  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max,
                                     levelindex_map<levelmultiindex, element_t> *lmi_map,
                                     levelindex_map<levelmultiindex, detail_t> &d_map) const
  {
    inverse_multiscale_transformation (l_min, l_max, lmi_map, d_map, mask);
  }

  /**
   * @brief Two-scale transform of one complete family (children -> parent + details)
   *
   * Pure computation, no map bookkeeping:
   *   u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
   *   d[k][i] = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
   *
   * @param data_on_siblings Element data of all NUM_CHILDREN children
   * @param data_on_coarse Output: parent element with u_coeffs, d_coeffs, vol, order
   * @param mask_coefficients Mask coefficient matrices M[k]
   */
  static void
  two_scale_family (const std::array<element_t, levelmultiindex::NUM_CHILDREN> &data_on_siblings,
                    detail_t &data_on_coarse, const std::vector<t8_mra::mat> &mask_coefficients)
  {
    const double scaling_factor = scaling_policy_t::forward_scaling_factor (levelmultiindex::NUM_CHILDREN);

    for (auto u = 0u; u < U_DIM; ++u) {
      std::array<double, DOF> u_parent;

      // Parent coefficients: u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
      for (auto i = 0u; i < DOF; ++i) {
        auto sum = 0.0;

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          const auto &Mk = mask_coefficients[k];
          const auto &uk = data_on_siblings[k].u_coeffs;
          for (auto j = 0u; j < DOF; ++j)
            sum += uk[element_t::dg_idx (u, j)] * Mk (j, i);
        }

        u_parent[i] = sum * scaling_factor;
        data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = u_parent[i];
      }

      // Detail coefficients: d[k][i] = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        const auto &Mk = mask_coefficients[k];
        const auto &uk = data_on_siblings[k].u_coeffs;
        for (auto i = 0u; i < DOF; ++i) {
          auto sum = 0.0;
          for (auto j = 0u; j < DOF; ++j)
            sum += Mk (i, j) * u_parent[j];

          data_on_coarse.d_coeffs[detail_t::wavelet_idx (k, u, i)] = uk[element_t::dg_idx (u, i)] - sum;
        }
      }
    }

    data_on_coarse.vol = data_on_siblings[0].vol * levelmultiindex::NUM_CHILDREN;
    data_on_coarse.order = data_on_siblings[0].order;
    ordering_policy_t::adjust_parent_order (data_on_coarse);
  }

  /** @brief Detail 2-norm per component; triangle scales by 1/sqrt(vol). */
  static std::array<double, U_DIM>
  detail_norm (const detail_t &detail)
  {
    std::array<double, U_DIM> norm = {};
    const auto &details = detail.d_coeffs;

    for (auto u = 0u; u < U_DIM; ++u) {
      double norm_sq = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
        for (auto i = 0u; i < DOF; ++i) {
          const auto d = details[detail_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }

      norm[u] = std::sqrt (norm_sq * scaling_policy_t::detail_norm_scale (detail.vol));
    }

    return norm;
  }

  /**
   * @brief Forward multiscale transformation (analysis), non-destructive.
   *
   * For every complete family in (l_min, l_max] computes the parent
   * (u_coeffs + details) via the two-scale relation and stores it in d_map;
   * lmi_map keeps its single-scale leaves. This is the decomposition step
   * coarsen and refine run to obtain the details their criteria act on. The
   * destructive full collapse is multiscale_decomposition.
   */
  static void
  multiscale_transformation (unsigned int l_min, unsigned int l_max,
                             levelindex_map<levelmultiindex, element_t> *lmi_map,
                             levelindex_map<levelmultiindex, detail_t> &d_map,
                             const std::vector<t8_mra::mat> &mask_coefficients)
  {
    index_set I_set;
    detail_t data_on_coarse;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_max; l > l_min; --l) {
      I_set.reserve (lmi_map->size (l));
      d_map[l - 1].reserve (lmi_map->size (l));

      for (const auto &[lmi, _] : lmi_map->operator[] (l))
        I_set.emplace (t8_mra::parent_lmi (lmi));

      for (const auto &lmi : I_set) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // Incomplete families (siblings on finer levels) carry no detail
        // information.
        auto family_complete = true;
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          const auto *sibling = lmi_map->find (siblings_lmi[k]);
          if (sibling == nullptr) {
            family_complete = false;
            break;
          }
          data_on_siblings[k] = *sibling;
        }
        if (!family_complete)
          continue;

        two_scale_family (data_on_siblings, data_on_coarse, mask_coefficients);
        d_map.insert (lmi, data_on_coarse);
      }

      I_set.clear ();
    }
  }

  /**
   * @brief Full forward multiscale transformation (restriction), destructive.
   *
   * Collapses the single-scale representation down to l_min: each complete
   * family is replaced by its parent in lmi_map (children erased) and its
   * details go to d_map. The single-scale leaves below l_min are gone; the
   * data now lives as (coarsest scaling coeffs + details). Inverted exactly by
   * inverse_multiscale_transformation. Used by the MST round-trip test; the
   * non-destructive analysis used during adaptation is multiscale_transformation.
   */
  static void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max, levelindex_map<levelmultiindex, element_t> *lmi_map,
                            levelindex_map<levelmultiindex, detail_t> &d_map,
                            const std::vector<t8_mra::mat> &mask_coefficients)
  {
    index_set I_set;
    detail_t data_on_coarse;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    for (auto l = l_max; l > l_min; --l) {
      // Collect all parent indices at level l-1
      for (const auto &[lmi, _] : lmi_map->operator[] (l))
        I_set.emplace (t8_mra::parent_lmi (lmi));

      d_map[l - 1].reserve (lmi_map->size (l));

      for (const auto &lmi : I_set) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // On an adaptive grid a family may be incomplete: some siblings stayed
        // refined on finer levels. Such families cannot be two-scale transformed;
        // their members remain leaves at level l.
        auto family_complete = true;
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          const auto *sibling = lmi_map->find (siblings_lmi[k]);
          if (sibling == nullptr) {
            family_complete = false;
            break;
          }
          data_on_siblings[k] = *sibling;
        }
        if (!family_complete)
          continue;

        two_scale_family (data_on_siblings, data_on_coarse, mask_coefficients);

        // The lmi_map leaf keeps only single-scale data (slice off d_coeffs).
        lmi_map->insert (lmi, static_cast<const element_t &> (data_on_coarse));
        d_map.insert (lmi, data_on_coarse);

        // Consume only this family's children; members of skipped (incomplete)
        // families must stay in the map as leaves.
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          lmi_map->erase (siblings_lmi[k]);
      }

      I_set.clear ();
    }
  }

  /**
   * @brief Inverse multiscale transformation (prolongation: coarse -> fine)
   *
   * Reconstructs children from parent and detail coefficients:
   *   u_child[k][i] = d[k][i] + Σ_j M[k](i,j) * u_parent[j]
   *
   * @param l_min Minimum refinement level
   * @param l_max Maximum refinement level
   * @param lmi_map Map from levelmultiindex to element data
   * @param d_map Detail coefficient storage
   * @param mask_coefficients Mask coefficient matrices M[k]
   */
  static void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max,
                                     levelindex_map<levelmultiindex, element_t> *lmi_map,
                                     levelindex_map<levelmultiindex, detail_t> &d_map,
                                     const std::vector<t8_mra::mat> &mask_coefficients)
  {
    element_t new_data;
    const double inv_scaling_factor = scaling_policy_t::inverse_scaling_factor ();

    for (auto l = l_min; l < l_max; ++l) {
      lmi_map->operator[] (l + 1).reserve (d_map[l].size ());

      for (const auto &[lmi, d] : d_map[l]) {
        const auto children_lmi = t8_mra::children_lmi (lmi);
        const auto lmi_data = lmi_map->get (lmi);
        const auto &details = d.d_coeffs;
        const auto &u_parent = lmi_data.u_coeffs;

        // Inverse MST: Reconstruct children u_child[k][i] = d[k][i] + Σ_j M[k](i,j) * u_parent[j]
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          const auto &Mk = mask_coefficients[k];
          for (auto u = 0u; u < U_DIM; ++u) {
            for (auto i = 0u; i < DOF; ++i) {
              auto sum = 0.0;

              for (auto j = 0u; j < DOF; ++j)
                sum += u_parent[element_t::dg_idx (u, j)] * Mk (i, j);

              new_data.u_coeffs[element_t::dg_idx (u, i)]
                = details[detail_t::wavelet_idx (k, u, i)] + sum * inv_scaling_factor;
            }
          }

          // Apply element-specific ordering adjustments
          new_data.vol = lmi_data.vol / levelmultiindex::NUM_CHILDREN;
          ordering_policy_t::adjust_child_order (new_data, k, lmi_data);

          lmi_map->insert (children_lmi[k], new_data);
        }

        lmi_map->erase (lmi);
      }

      d_map.erase (l);
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
