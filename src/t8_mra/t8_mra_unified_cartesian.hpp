#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_base.hpp"
#include "t8_mra/t8_mra_adaptation.hpp"
#include "t8_mra/num/mask_coefficients_compute.hxx"

namespace t8_mra
{

/**
 * @brief Cartesian multiscale implementation (partial specialization)
 *
 * Generic implementation for all cartesian elements (LINE, QUAD, HEX).
 * Inherits all common MRA functionality from multiscale_base and
 * implements cartesian-specific:
 *   - Detail norm computation (uses sqrt(volume) scaling)
 *   - Projection (uses Gauss-Legendre tensor product quadrature)
 *   - No vertex ordering needed (cartesian elements have standard ordering)
 *
 * @tparam TShape Element shape (T8_ECLASS_LINE, T8_ECLASS_QUAD, or T8_ECLASS_HEX)
 * @tparam U Number of solution components
 * @tparam P Polynomial order
 */
template <t8_eclass TShape, int U, int P>
  requires is_cartesian<TShape>
class multiscale<TShape, U, P>:
  public multiscale_base<TShape, U, P>,
  public multiscale_adaptation<multiscale<TShape, U, P>> {

 public:
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
   * @brief Constructor for cartesian multiscale
   *
   * @param _max_level Maximum refinement level
   * @param _c_thresh Coarsening threshold
   * @param _gamma Gamma parameter for thresholding
   * @param _num_quad_points_1d Number of 1D Gauss-Legendre quadrature points
   * @param _balanced Whether to balance the forest
   * @param _comm MPI communicator
   */
  multiscale (int _max_level, double _c_thresh, int _gamma, int _num_quad_points_1d, bool _balanced, sc_MPI_Comm _comm)
    : Base (_max_level, _c_thresh, _gamma, _num_quad_points_1d, _balanced, _comm)
  {
    // Initialize mask coefficients via computation
    t8_mra::initialize_mask_coefficients_computed<TShape> (Base::P_DIM, Base::DOF, Base::mask_coefficients,
                                                           Base::inverse_mask_coefficients);
  }

  //=============================================================================
  // Element-Specific Implementation: Detail Norm
  //=============================================================================

  /**
   * @brief Compute local detail norm for cartesian elements
   *
   * Cartesian elements use sqrt(volume) scaling factor (no factor of 2).
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

    // Cartesian-specific scaling: sqrt(volume) without factor of 2
    const auto scaling_factor = std::sqrt (vol);

    for (auto u = 0u; u < Base::U_DIM; ++u) {
      double norm_sq = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        for (auto i = 0u; i < Base::DOF; ++i) {
          const auto d = details[element_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }
      }
      detail_norm[u] = std::sqrt (norm_sq) / scaling_factor;
    }

    return detail_norm;
  }

  //=============================================================================
  // Element-Specific Implementation: Projection
  //=============================================================================

  /**
   * @brief Project a function onto the DG basis for a cartesian element
   *
   * Uses Gauss-Legendre tensor product quadrature on reference element [0,1]^DIM.
   * Direct mapping from physical to reference coordinates (no transformation matrix).
   *
   * @param dg_coeffs Output vector for DG coefficients
   * @param tree_idx Tree index in forest
   * @param element Pointer to element
   * @param func Function to project
   */
  template <typename Func>
  void
  project_impl (std::vector<double> &dg_coeffs, int tree_idx, const t8_element_t *element, Func &&func)
  {
    // Extract element vertices (cartesian: axis-aligned box)
    constexpr int num_vertices = (Base::DIM == 1 ? 2 : (Base::DIM == 2 ? 4 : 8));
    double vertices[8][3] = {};

    if constexpr (Base::DIM == 2 && TShape == T8_ECLASS_QUAD) {
      // t8code QUAD vertex ordering: 0-1-2-3 as (0,0)-(1,0)-(0,1)-(1,1)
      // We need: 0-1-2-3 as (0,0)-(1,0)-(1,1)-(0,1) for standard quad
      // So we swap vertices 2 and 3
      const int vertex_perm[4] = { 0, 1, 3, 2 };
      for (int i = 0; i < num_vertices; ++i)
        t8_forest_element_coordinate (Base::forest, tree_idx, element, vertex_perm[i], vertices[i]);
    }
    else {
      for (int i = 0; i < num_vertices; ++i)
        t8_forest_element_coordinate (Base::forest, tree_idx, element, i, vertices[i]);
    }

    // Get physical quadrature points via direct mapping
    const auto phys_quad_points = Base::DG_basis.deref_quad_points (vertices);

    // ONE-TO-ONE implementation from old t8_mra_cartesian.hpp::project()
    // Note: For orthonormal Legendre basis, NO volume/Jacobian scaling in projection!
    // The basis functions are normalized on reference element [0,1]^DIM

    // Precompute basis values at all quadrature points
    std::vector<std::array<double, Base::DOF>> basis_at_quad (Base::DG_basis.num_quad_points);
    for (auto q = 0u; q < Base::DG_basis.num_quad_points; ++q) {
      // Extract reference quadrature point
      std::vector<double> x_ref (Base::DIM);
      for (unsigned int d = 0; d < Base::DIM; ++d)
        x_ref[d] = Base::DG_basis.ref_quad_points[Base::DIM * q + d];

      const auto basis_vals = Base::DG_basis.basis_value (x_ref);
      for (auto i = 0u; i < Base::DOF; ++i)
        basis_at_quad[q][i] = basis_vals[i];
    }

    // Project onto each basis function
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto q = 0u; q < Base::DG_basis.num_quad_points; ++q) {
        // Extract physical coordinates
        std::array<double, Base::DIM> x_phys;
        for (unsigned int d = 0; d < Base::DIM; ++d)
          x_phys[d] = phys_quad_points[Base::DIM * q + d];

        // Evaluate function at physical point
        // Supports two calling conventions:
        //   func(x, y, z) -> std::array<double, U>  (return value)
        //   func(x, y, z, double* out)               (output pointer)
        std::array<double, Base::U_DIM> f_val;
        if constexpr (Base::DIM == 1) {
          if constexpr (std::is_invocable_v<decltype (func), double>)
            f_val = func (x_phys[0]);
          else
            func (x_phys[0], f_val.data ());
        }
        else if constexpr (Base::DIM == 2) {
          if constexpr (std::is_invocable_v<decltype (func), double, double>)
            f_val = func (x_phys[0], x_phys[1]);
          else
            func (x_phys[0], x_phys[1], f_val.data ());
        }
        else {
          if constexpr (std::is_invocable_v<decltype (func), double, double, double>)
            f_val = func (x_phys[0], x_phys[1], x_phys[2]);
          else
            func (x_phys[0], x_phys[1], x_phys[2], f_val.data ());
        }

        // Accumulate quadrature sum: integral(f * phi_i)
        // Note: For orthonormal basis, the volume scaling cancels out
        for (auto u = 0u; u < Base::U_DIM; ++u)
          sum[u] += Base::DG_basis.quad_weights[q] * f_val[u] * basis_at_quad[q][i];
      }

      // Store coefficients directly (already includes volume scaling from quadrature)
      for (auto u = 0u; u < Base::U_DIM; ++u)
        dg_coeffs[element_t::dg_idx (u, i)] = sum[u];
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
  template <typename Func>
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, Func &&func)
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

        data_element.vol = t8_forest_element_volume (Base::forest, tree_idx, element);

        // Project function using element-specific projection
        project_impl (data_element.u_coeffs, tree_idx, element, func);

        user_data->lmi_map->insert (lmi, data_element);

        // Insert lmi into forest
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
      }
    }

    T8_FREE (elem_data);
    t8_forest_set_user_data (Base::forest, user_data);
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
