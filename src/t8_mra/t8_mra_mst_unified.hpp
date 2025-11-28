#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/data/levelindex_map.hpp"
#include "t8_mra/data/levelindex_set.hpp"
#include "t8_mra/data/triangle_order.hpp"
#include "t8_mra/num/mat.hpp"
#include "t8_mra/t8_basis.hpp"

#include <array>
#include <iostream>

namespace t8_mra
{

/**
 * @brief Ordering policy for element-specific vertex handling
 *
 * Default policy: no-op for cartesian elements (QUAD, LINE, HEX)
 */
template <t8_eclass TShape>
struct OrderingPolicy
{
  /**
   * @brief Adjust parent element ordering (no-op for cartesian elements)
   */
  template <typename T>
  static void
  adjust_parent_order (T &data)
  {
    // No ordering adjustment needed for cartesian elements
  }

  /**
   * @brief Adjust child element ordering (no-op for cartesian elements)
   */
  template <typename T>
  static void
  adjust_child_order (T &child_data, int child_id, const T &parent_data)
  {
    // No ordering adjustment needed for cartesian elements
  }
};

/**
 * @brief Ordering policy specialization for triangles
 *
 * Triangles require complex vertex ordering to maintain consistency
 * across refinement levels
 */
template <>
struct OrderingPolicy<T8_ECLASS_TRIANGLE>
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
struct MSTScalingPolicy
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
struct MSTScalingPolicy<T8_ECLASS_TRIANGLE>
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
 * @brief Unified multiscale transformation (MST) operations
 *
 * This class template provides unified forward and inverse multiscale
 * transformations that work for both triangular and cartesian elements.
 * Element-specific behavior is controlled via policy classes.
 *
 * @tparam TElement Element type (data_per_element<TShape, U, P>)
 * @tparam OrderingPolicyT Policy for vertex ordering (default: OrderingPolicy<TElement::Shape>)
 * @tparam ScalingPolicyT Policy for scaling factors (default: MSTScalingPolicy<TElement::Shape>)
 */
template <typename TElement, typename OrderingPolicyT = OrderingPolicy<TElement::Shape>,
          typename ScalingPolicyT = MSTScalingPolicy<TElement::Shape>>
class UnifiedMST
{
 public:
  using element_t = TElement;
  using levelmultiindex = t8_mra::levelmultiindex<TElement::Shape>;
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;

  static constexpr auto Shape = TElement::Shape;
  static constexpr unsigned int U_DIM = TElement::U_DIM;
  static constexpr unsigned int DOF = TElement::DOF;

  /**
   * @brief Forward multiscale transformation (restriction: fine -> coarse)
   *
   * Computes parent coefficients and detail coefficients:
   *   u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
   *   d[k][i] = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
   *
   * @param l_min Minimum refinement level
   * @param l_max Maximum refinement level
   * @param lmi_map Map from levelmultiindex to element data
   * @param d_map Detail coefficient storage
   * @param mask_coefficients Mask coefficient matrices M[k]
   */
  static void
  forward_transformation (unsigned int l_min, unsigned int l_max,
                          levelindex_map<levelmultiindex, element_t> *lmi_map,
                          levelindex_map<levelmultiindex, element_t> &d_map,
                          const std::vector<t8_mra::mat> &mask_coefficients)
  {
    index_set I_set;
    element_t data_on_coarse;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    const double scaling_factor = ScalingPolicyT::forward_scaling_factor (levelmultiindex::NUM_CHILDREN);

    for (auto l = l_max; l > l_min; --l) {
      // Collect all parent indices at level l-1
      for (const auto &[lmi, _] : lmi_map->operator[] (l))
        I_set.emplace (t8_mra::parent_lmi (lmi));

      d_map[l - 1].reserve (lmi_map->size (l));

      for (const auto &lmi : I_set) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // Load children - LMI structure encodes the ordering
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          data_on_siblings[k] = lmi_map->get (siblings_lmi[k]);

        // Forward MST: Compute parent coefficients
        // u_parent[i] = scaling * Σ_k Σ_j u_child[k][j] * M[k](j,i)
        for (auto u = 0u; u < U_DIM; ++u) {
          for (auto i = 0u; i < DOF; ++i) {
            auto sum = 0.0;

            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
              for (auto j = 0u; j < DOF; ++j)
                sum += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i);

            data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = sum * scaling_factor;
            data_on_coarse.vol = data_on_siblings[0].vol * levelmultiindex::NUM_CHILDREN;
          }

          // Compute detail coefficients: d[k][i] = u_child[k][i] - Σ_j M[k](i,j) * u_parent[j]
          for (auto i = 0u; i < DOF; ++i) {
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
              auto sum = 0.0;
              for (auto j = 0u; j < DOF; ++j)
                sum += mask_coefficients[k](i, j) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];

              data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, i)]
                = data_on_siblings[k].u_coeffs[element_t::dg_idx (u, i)] - sum;
            }
          }
        }

        // Apply element-specific ordering adjustments
        data_on_coarse.order = data_on_siblings[0].order;
        OrderingPolicyT::adjust_parent_order (data_on_coarse);

        lmi_map->insert (lmi, data_on_coarse);
        d_map[l - 1].emplace (lmi, data_on_coarse);
      }

      lmi_map->erase (l);
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
  inverse_transformation (unsigned int l_min, unsigned int l_max,
                          levelindex_map<levelmultiindex, element_t> *lmi_map,
                          levelindex_map<levelmultiindex, element_t> &d_map,
                          const std::vector<t8_mra::mat> &mask_coefficients)
  {
    element_t new_data;
    const double inv_scaling_factor = ScalingPolicyT::inverse_scaling_factor ();

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
                = details[element_t::wavelet_idx (k, u, i)] + sum * inv_scaling_factor;
            }
          }

          // Apply element-specific ordering adjustments
          new_data.vol = lmi_data.vol / levelmultiindex::NUM_CHILDREN;
          OrderingPolicyT::adjust_child_order (new_data, k, lmi_data);

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
