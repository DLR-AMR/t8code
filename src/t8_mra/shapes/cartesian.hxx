#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/base.hxx"
#include "t8_mra/core/adaptation.hxx"
#include "t8_mra/num/mask_coefficients_compute.hxx"

#include <span>

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
  using Adaptation::coarsen;
  using Adaptation::refine;
  using Adaptation::initialize_data_adaptive;
  using Adaptation::balance;

  //=============================================================================
  // Constructor
  //=============================================================================

  /**
   * @brief Constructor for cartesian multiscale
   *
   * @param _max_level Maximum refinement level
   * @param _comm MPI communicator
   */
  multiscale (int _max_level, sc_MPI_Comm _comm): Base (_max_level, _comm)
  {
    // Initialize mask coefficients via computation
    t8_mra::initialize_mask_coefficients_computed<TShape> (Base::P_DIM, Base::DOF, Base::mask_coefficients);
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
          const auto d = details[Base::detail_t::wavelet_idx (k, u, i)];
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
  project_impl (std::span<double> dg_coeffs, int tree_idx, const t8_element_t *element, Func &&func)
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
    const auto phys_quad_points = Base::basis.deref_quad_points (vertices);

    
    // Note: For orthonormal Legendre basis, NO volume/Jacobian scaling in projection!
    // The basis functions are normalized on reference element [0,1]^DIM

    // Precompute basis values at all quadrature points
    std::vector<std::array<double, Base::DOF>> basis_at_quad (Base::basis.quad.num_points);
    for (auto q = 0u; q < Base::basis.quad.num_points; ++q) {
      // Extract reference quadrature point
      std::vector<double> x_ref (Base::DIM);
      for (unsigned int d = 0; d < Base::DIM; ++d)
        x_ref[d] = Base::basis.quad.points[Base::DIM * q + d];

      const auto basis_vals = Base::basis.basis_value (x_ref);
      for (auto i = 0u; i < Base::DOF; ++i)
        basis_at_quad[q][i] = basis_vals[i];
    }

    // Project onto each basis function
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto q = 0u; q < Base::basis.quad.num_points; ++q) {
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
          sum[u] += Base::basis.quad.weights[q] * f_val[u] * basis_at_quad[q][i];
      }

      // Store coefficients directly (already includes volume scaling from quadrature)
      for (auto u = 0u; u < Base::U_DIM; ++u)
        dg_coeffs[element_t::dg_idx (u, i)] = sum[u];
    }
  }

  /**
   * @brief Project a function onto a single forest leaf
   *
   * Builds the complete element data (volume, DG coefficients) for the
   * given leaf.
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

    data.vol = t8_forest_element_volume (Base::forest, tree_idx, element);
    project_impl (data.u_coeffs, tree_idx, element, func);

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
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
