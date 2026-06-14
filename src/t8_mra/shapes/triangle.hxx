#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/base.hxx"
#include "t8_mra/core/adaptation.hxx"
#include "t8_mra/data/triangle_order.hxx"
#include "t8_mra/num/mask_coefficients.hxx"

#include <span>

namespace t8_mra
{

/**
 * @brief Triangle-specific multiscale implementation (partial specialization)
 *
 * Inherits all common MRA functionality from multiscale_base and
 * implements triangle-specific:
 *   - Detail norm computation (uses sqrt(2*volume) scaling)
 *   - Projection (uses Dunavant quadrature + transformation matrix)
 *   - Vertex ordering (via triangle_order)
 */
template <int U, int P>
class multiscale<T8_ECLASS_TRIANGLE, U, P>:
  public multiscale_base<T8_ECLASS_TRIANGLE, U, P>,
  public multiscale_adaptation<multiscale<T8_ECLASS_TRIANGLE, U, P>> {

 public:
  static constexpr t8_eclass TShape = T8_ECLASS_TRIANGLE;
  using Base = multiscale_base<TShape, U, P>;
  using Adaptation = multiscale_adaptation<multiscale<TShape, U, P>>;
  using element_t = typename Base::element_t;
  using levelmultiindex = typename Base::levelmultiindex;
  using index_set = typename Base::index_set;

  // Make adaptation methods accessible
  using Adaptation::coarsen;
  using Adaptation::refine;
  using Adaptation::initialize_data_adaptive;
  using Adaptation::balance;

  //=============================================================================
  // Constructor
  //=============================================================================

  /**
   * @brief Constructor for triangle multiscale
   *
   * @param _max_level Maximum refinement level
   * @param _comm MPI communicator
   */
  multiscale (int _max_level, sc_MPI_Comm _comm): Base (_max_level, _comm)
  {
    // Initialize mask coefficients from hardcoded values
    t8_mra::initialize_mask_coefficients<TShape> (Base::P_DIM, Base::DOF, Base::mask_coefficients,
                                                  Base::inverse_mask_coefficients);
  }

  //=============================================================================
  // Element-Specific Implementation: Detail Norm
  //=============================================================================

  /**
   * @brief Compute local detail norm for triangles
   *
   * Triangles use sqrt(2*volume) scaling factor in the detail norm computation.
   *
   * @param lmi Level multi-index
   * @return Array of detail norms (one per solution component)
   */
  std::array<double, Base::U_DIM>
  local_detail_norm (const levelmultiindex &lmi) override
  {
    std::array<double, Base::U_DIM> detail_norm = {};
    const auto vol = Base::d_map.get (lmi).vol;
    const auto &details = Base::d_map.get (lmi).d_coeffs;

    for (auto u = 0u; u < Base::U_DIM; ++u) {
      double norm_sq = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        for (auto i = 0u; i < Base::DOF; ++i) {
          const auto d = details[Base::detail_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }
      }
      // Match old implementation: sqrt(sum / vol)
      detail_norm[u] = std::sqrt (norm_sq / vol);
    }

    return detail_norm;
  }

  //=============================================================================
  // Element-Specific Implementation: Projection
  //=============================================================================

  /**
   * @brief Project a function onto the DG basis for a triangle
   *
   * Uses Dunavant quadrature and requires transformation to reference element
   * due to arbitrary triangle geometry.
   *
   * @param dg_coeffs Output vector for DG coefficients
   * @param tree_idx Tree index in forest
   * @param element Pointer to element
   * @param func Function to project (signature: std::array<double, U_DIM>(double x, double y))
   */
  template <typename Func>
  void
  project_impl (std::span<double> dg_coeffs, int tree_idx, const t8_element_t *element,
                const std::array<int, 3> &point_order, Func &&func)
  {
    // Reorder t8code vertices into reference order: t8code vertex i goes to
    // position point_order[i].
    double vertices[3][3];
    for (auto i = 0; i < 3; ++i)
      t8_forest_element_coordinate (Base::forest, tree_idx, element, i, vertices[point_order[i]]);

    // Compute transformation to reference element
    auto [trafo_mat, perm] = Base::basis.trafo_matrix_to_ref_element (vertices);
    const auto deref_quad_points = Base::basis.deref_quad_points (vertices);
    const auto volume = t8_forest_element_volume (Base::forest, tree_idx, element);

    const auto scaling_factor = shape_traits<element_t::Shape>::basis_normalization (volume);

    // Project function onto DG basis using Dunavant quadrature
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto j = 0u; j < Base::basis.num_quad_points; ++j) {
        const auto x_deref = deref_quad_points[2 * j];
        const auto y_deref = deref_quad_points[1 + 2 * j];

        // Transform to reference coordinates
        const auto ref = Base::basis.ref_point (trafo_mat, perm, { x_deref, y_deref, 1.0 });

        // Evaluate function at physical point
        const auto f_val = func (x_deref, y_deref);

        // Evaluate basis at reference point
        const auto basis_val = Base::basis.basis_value (ref);

        // Accumulate: ∫ f(x) φᵢ(x) dx ≈ Σ w_j f(x_j) φᵢ(x_j) * scaling
        for (auto k = 0u; k < Base::U_DIM; ++k)
          sum[k] += Base::basis.quad_weights[j] * f_val[k] * scaling_factor * basis_val[i];
      }

      for (auto k = 0u; k < Base::U_DIM; ++k)
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k] * volume;
    }
  }

  /**
   * @brief Project a function onto a single forest leaf
   *
   * Builds the complete element data (vertex order, volume, DG coefficients)
   * for the given leaf.
   *
   * @param tree_idx Tree index in forest
   * @param element Pointer to element
   * @param func Function to project
   * @return Element data of the leaf
   */
  template <typename Func>
  element_t
  project_leaf (int tree_idx, const t8_element_t *element, Func &&func)
  {
    element_t data;

    const auto gtreeid = t8_forest_global_tree_id (Base::forest, tree_idx);
    const auto *scheme = t8_forest_get_scheme (Base::forest);
    triangle_order::get_point_order_at_level (gtreeid, element, scheme, data.order);

    data.vol = t8_forest_element_volume (Base::forest, tree_idx, element);
    project_impl (data.u_coeffs, tree_idx, element, data.order, func);

    return data;
  }

  //=============================================================================
  // Initialization
  //=============================================================================

  /**
   * @brief Initialize data on uniform mesh with function projection
   *
   * @param mesh Coarse mesh
   * @param scheme Element scheme
   * @param level Uniform refinement level
   * @param func Function to project onto the mesh
   */
  template <typename Func>
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, Func &&func)
  {
    Base::forest = t8_forest_new_uniform (mesh, scheme, level, 0, Base::comm);
    Base::build_lmi_map (scheme, [&] (int tree_idx, const t8_element_t *element) {
      return project_leaf (tree_idx, element, func);
    });
  }

  //=============================================================================
  // Post-Adaptation Hook
  //=============================================================================

  /**
   * @brief Update vertex orders after forest adaptation
   *
   * Triangle vertex order depends on ancestor.type which is not stored in LMI.
   * This function updates the vertex orders in lmi_map to match the actual
   * t8code elements after adaptation.
   */
  void
  post_adaptation_hook () override
  {
    if (Base::forest == nullptr)
      return;

    auto *forest_data = Base::get_user_data ();
    const auto *forest_scheme = t8_forest_get_scheme (Base::forest);

    Base::for_each_local_leaf (
      [&] (t8_locidx_t tree_idx, const t8_element_t *elem, unsigned int local_idx, t8_gloidx_t) {
        const auto lmi = t8_mra::get_lmi_from_forest_data<element_t> (forest_data, local_idx);
        if (auto *elem_data = forest_data->lmi_map->find (lmi))
          t8_mra::triangle_order::get_point_order_at_level (tree_idx, elem, forest_scheme, elem_data->order);
      });
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
