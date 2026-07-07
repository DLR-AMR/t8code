#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelindex_set.hxx"
#include "t8_mra/data/triangle_order.hxx"
#include "t8_mra/num/mat.hxx"
#include "t8_mra/num/dg_basis.hxx"

#include <algorithm>
#include <array>

namespace t8_mra
{

/**
 * @brief Ordering policy for element-specific vertex handling
 *
 * Default policy: no-op for cartesian elements (QUAD, LINE, HEX)
 */
template <t8_eclass TShape>
struct ordering_policy
{
  /**
   * @brief Adjust parent element ordering (no-op for cartesian elements)
   */
  template <typename T>
  static void
  adjust_parent_order (T &data)
  {
  }

  /**
   * @brief Adjust child element ordering (no-op for cartesian elements)
   */
  template <typename T>
  static void
  adjust_child_order (T &child_data, int child_id, const T &parent_data)
  {
  }
};

/**
 * @brief Ordering policy specialization for triangles
 *
 * Triangles require complex vertex ordering to maintain consistency
 * across refinement levels
 */
template <>
struct ordering_policy<T8_ECLASS_TRIANGLE>
{
  /**
   * @brief Compute parent vertex order from child order
   */
  template <typename T>
  static void
  adjust_parent_order (T &data)
  {
    triangle_order::get_parent_order (data.order);
  }

  /**
   * @brief Compute child vertex order from parent order
   */
  template <typename T>
  static void
  adjust_child_order (T &child_data, int child_id, const T &parent_data)
  {
    child_data.order = parent_data.order;
    triangle_order::get_point_order (child_data.order, child_id);
  }
};

/**
 * @brief Scaling policy for MST normalization
 *
 * Different element types may require different scaling factors
 * in the multiscale transformation
 */
template <t8_eclass TShape>
struct mst_scaling_policy
{
  /**
   * @brief Forward MST normalization factor
   *
   * For cartesian elements with L² orthonormal basis on reference element,
   * we need to average over children: factor = 1/NUM_CHILDREN
   */
  static constexpr double
  forward_scaling_factor (unsigned int num_children)
  {
    return 1.0 / static_cast<double> (num_children);
  }

  /**
   * @brief Inverse MST scaling factor (typically 1.0)
   */
  static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }
};

/**
 * @brief Scaling policy specialization for triangles
 *
 * Triangles use a different normalization convention
 */
template <>
struct mst_scaling_policy<T8_ECLASS_TRIANGLE>
{
  /**
   * @brief Forward MST normalization factor for triangles
   *
   * Triangle implementation does NOT divide by NUM_CHILDREN
   */
  static constexpr double
  forward_scaling_factor (unsigned int /* num_children */)
  {
    return 1.0;
  }

  /**
   * @brief Inverse MST scaling factor
   */
  static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }
};

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

        // Incomplete families (siblings on finer levels) carry no detail information
        const auto family_complete
          = std::all_of (siblings_lmi.begin (), siblings_lmi.end (),
                         [&] (const levelmultiindex &sibling) { return lmi_map->contains (sibling); });
        if (!family_complete)
          continue;

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          data_on_siblings[k] = lmi_map->get (siblings_lmi[k]);

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
        const auto family_complete
          = std::all_of (siblings_lmi.begin (), siblings_lmi.end (),
                         [&] (const levelmultiindex &sibling) { return lmi_map->contains (sibling); });
        if (!family_complete)
          continue;

        // Load children - LMI structure encodes the ordering
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          data_on_siblings[k] = lmi_map->get (siblings_lmi[k]);

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

        // Inverse MST: Reconstruct children u_child[k][i] = d[k][i] + Σ_j M[k](i,j) * u_parent[j]
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          for (auto u = 0u; u < U_DIM; ++u) {
            for (auto i = 0u; i < DOF; ++i) {
              auto sum = 0.0;

              for (auto j = 0u; j < DOF; ++j)
                sum += lmi_data.u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j);

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
