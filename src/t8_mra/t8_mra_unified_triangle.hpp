#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_base.hpp"
#include "t8_mra/t8_mra_adaptation.hpp"
#include "t8_mra/data/triangle_order.hpp"
#include "t8_mra/num/mask_coefficients.hpp"

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
  using Adaptation::coarsening_new;
  using Adaptation::refinement_new;

  //=============================================================================
  // Constructor
  //=============================================================================

  /**
   * @brief Constructor for triangle multiscale
   *
   * @param _max_level Maximum refinement level
   * @param _c_thresh Coarsening threshold
   * @param _gamma Gamma parameter for thresholding
   * @param _dunavant_rule Dunavant quadrature rule
   * @param _balanced Whether to balance the forest
   * @param _comm MPI communicator
   */
  multiscale (int _max_level, double _c_thresh, int _gamma, int _dunavant_rule, bool _balanced, sc_MPI_Comm _comm)
    : Base (_max_level, _c_thresh, _gamma, _dunavant_rule, _balanced, _comm)
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
          const auto d = details[element_t::wavelet_idx (k, u, i)];
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
  project_impl (std::vector<double> &dg_coeffs, int tree_idx, const t8_element_t *element,
                const std::array<int, 3> &point_order, Func &&func)
  {
    // ONE-TO-ONE implementation of old t8_mra.hpp::project()
    // Takes point_order as parameter (already computed in initialize_data)

    // Get triangle vertices from t8code and reorder them
    // OLD IMPLEMENTATION: vertices[order[i]] = t8code_vertex[i]
    // This stores t8code vertex i at position order[i] in the array
    double vertices[3][3];
    for (auto i = 0; i < 3; ++i)
      t8_forest_element_coordinate (Base::forest, tree_idx, element, i, vertices[point_order[i]]);

    // Compute transformation to reference element
    auto [trafo_mat, perm] = Base::DG_basis.trafo_matrix_to_ref_element (vertices);
    const auto deref_quad_points = Base::DG_basis.deref_quad_points (vertices);
    const auto volume = t8_forest_element_volume (Base::forest, tree_idx, element);

    // Triangle-specific scaling factor
    const auto scaling_factor = std::sqrt (1.0 / (2.0 * volume));

    // Project function onto DG basis using Dunavant quadrature
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto j = 0u; j < Base::DG_basis.num_quad_points; ++j) {
        const auto x_deref = deref_quad_points[2 * j];
        const auto y_deref = deref_quad_points[1 + 2 * j];

        // Transform to reference coordinates
        const auto ref = Base::DG_basis.ref_point (trafo_mat, perm, { x_deref, y_deref, 1.0 });

        // Evaluate function at physical point
        const auto f_val = func (x_deref, y_deref);

        // Evaluate basis at reference point
        const auto basis_val = Base::DG_basis.basis_value (ref);

        // Accumulate: ∫ f(x) φᵢ(x) dx ≈ Σ w_j f(x_j) φᵢ(x_j) * scaling
        for (auto k = 0u; k < Base::U_DIM; ++k)
          sum[k] += Base::DG_basis.quad_weights[j] * f_val[k] * scaling_factor * basis_val[i];
      }

      for (auto k = 0u; k < Base::U_DIM; ++k)
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k] * volume;
    }
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
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, auto &&func)
  {
    Base::forest = t8_forest_new_uniform (mesh, scheme, level, 0, Base::comm);

    levelmultiindex *elem_data;
    t8_mra::forest_data<element_t> *user_data;

    user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    elem_data = T8_ALLOC (levelmultiindex, 1);

    T8_ASSERT (t8_forest_is_committed (Base::forest));

    const auto num_local_elements = t8_forest_get_global_num_leaf_elements (Base::forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (Base::forest);

    user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (Base::maximum_level);
    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);
    user_data->mra_instance = this;
    user_data->current_refinement_level = level;

    const auto num_local_trees = t8_forest_get_num_local_trees (Base::forest);
    auto current_idx = 0u;

    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (Base::forest, tree_idx);
      const auto base_element = t8_forest_global_tree_id (Base::forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements_in_tree; ++ele_idx, ++current_idx) {
        element_t data_element;
        const auto *element = t8_forest_get_leaf_element_in_tree (Base::forest, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (base_element, element, scheme);

        // Get vertex ordering for projection
        std::array<int, 3> point_order;
        triangle_order::get_point_order_at_level (base_element, element, scheme, point_order);

        data_element.order = point_order;
        data_element.vol = t8_forest_element_volume (Base::forest, tree_idx, element);

        // Project function using element-specific projection (pass point_order like old version)
        project_impl (data_element.u_coeffs, tree_idx, element, point_order, func);

        user_data->lmi_map->insert (lmi, data_element);

        // Insert lmi into forest
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
      }
    }

    T8_FREE (elem_data);
    t8_forest_set_user_data (Base::forest, user_data);
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
    t8_locidx_t elem_offset = 0;

    for (t8_locidx_t tree_idx = 0; tree_idx < t8_forest_get_num_local_trees (Base::forest); ++tree_idx) {
      const auto num_elems = t8_forest_get_tree_num_leaf_elements (Base::forest, tree_idx);
      const auto base_element = tree_idx;

      for (t8_locidx_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
        const auto *elem = t8_forest_get_leaf_element_in_tree (Base::forest, tree_idx, elem_idx);
        const auto lmi = t8_mra::get_lmi_from_forest_data<element_t> (forest_data, elem_offset + elem_idx);

        if (forest_data->lmi_map->contains (lmi)) {
          std::array<int, 3> correct_order;
          t8_mra::triangle_order::get_point_order_at_level (base_element, elem, forest_scheme, correct_order);

          auto &elem_data = forest_data->lmi_map->get (lmi);
          elem_data.order = correct_order;
        }
      }
      elem_offset += num_elems;
    }
  }

  //=============================================================================
  // Cleanup
  //=============================================================================

  /**
   * @brief Clean up forest data structures
   */
  void
  cleanup ()
  {
    Base::cleanup ();

    auto *user_data = Base::get_user_data ();
    if (user_data) {
      if (user_data->lmi_map) {
        delete user_data->lmi_map;
        user_data->lmi_map = nullptr;
      }
      if (user_data->lmi_idx) {
        sc_array_destroy (user_data->lmi_idx);
        user_data->lmi_idx = nullptr;
      }
      T8_FREE (user_data);
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
