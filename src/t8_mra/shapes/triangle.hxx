#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/base.hxx"
#include "t8_mra/core/adaptation.hxx"
#include "t8_mra/data/triangle_order.hxx"
#include "t8_mra/num/mask_coefficients.hxx"

#include <span>
#include <array>
#include <vector>

namespace t8_mra
{

/**
 * @brief Triangle-specific multiscale implementation (partial specialization)
 *
 * Inherits all common MRA functionality from multiscale_base and
 * implements triangle-specific:
 *   - Detail norm computation (1/sqrt(volume) scaling)
 *   - Projection (uses Dunavant quadrature + transformation matrix)
 *   - Vertex ordering (via triangle_order)
 */
template <int U, int P>
class multiscale<T8_ECLASS_TRIANGLE, U, P>:
  public multiscale_base<multiscale<T8_ECLASS_TRIANGLE, U, P>, T8_ECLASS_TRIANGLE, U, P>,
  public multiscale_adaptation<multiscale<T8_ECLASS_TRIANGLE, U, P>> {

 public:
  static constexpr t8_eclass TShape = T8_ECLASS_TRIANGLE;
  using Base = multiscale_base<multiscale<TShape, U, P>, TShape, U, P>;
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
    // Compute two-scale mask coefficients (reference basis + quadrature).
    t8_mra::compute_mask<TShape, Base::P_DIM> (Base::mask_coefficients);
  }

  //=============================================================================
  // Element-Specific Implementation: Detail Norm
  //=============================================================================

  /**
   * @brief Compute local detail norm for triangles
   *
   * Detail 2-norm over children scaled by 1/sqrt(volume).
   *
   * @param lmi Level multi-index
   * @return Array of detail norms (one per solution component)
   */
  std::array<double, Base::U_DIM>
  local_detail_norm (const levelmultiindex &lmi)
  {
    std::array<double, Base::U_DIM> detail_norm = {};
    const auto &data = Base::d_map.get (lmi);
    const auto vol = data.vol;
    const auto &details = data.d_coeffs;

    for (auto u = 0u; u < Base::U_DIM; ++u) {
      double norm_sq = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        for (auto i = 0u; i < Base::DOF; ++i) {
          const auto d = details[Base::detail_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }
      }
      detail_norm[u] = std::sqrt (norm_sq / vol);
    }

    return detail_norm;
  }

  //=============================================================================
  // Element-Specific Implementation: Projection
  //=============================================================================

  /// Corner coordinates in reference order: t8code vertex i goes to position
  /// order[i].
  void
  element_vertex_coords (int tree_idx, const t8_element_t *element, const std::array<int, 3> &order,
                         double vertices[3][3])
  {
    for (auto i = 0; i < 3; ++i)
      t8_forest_element_coordinate (Base::forest, tree_idx, element, i, vertices[order[i]]);
  }

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
    double vertices[3][3];
    element_vertex_coords (tree_idx, element, point_order, vertices);

    auto [trafo_mat, perm] = Base::basis.trafo_matrix_to_ref_element (vertices);
    const auto deref_quad_points = Base::basis.deref_quad_points (vertices);
    const auto volume = t8_forest_element_volume (Base::forest, tree_idx, element);

    const auto scaling_factor = basis<element_t::Shape, Base::P_DIM>::normalization (volume);

    // Precompute per quad point (independent of basis index i): all DOF basis
    // values and the function value.
    const auto num_q = Base::basis.quad.num_points;
    std::vector<std::array<double, Base::DOF>> basis_at_quad (num_q);
    std::vector<std::array<double, Base::U_DIM>> f_at_quad (num_q);
    for (auto j = 0u; j < num_q; ++j) {
      const auto x_deref = deref_quad_points[2 * j];
      const auto y_deref = deref_quad_points[1 + 2 * j];

      const auto ref = Base::basis.ref_point (trafo_mat, perm, { x_deref, y_deref });
      basis_at_quad[j] = Base::basis.basis_value ({ ref[0], ref[1] });
      f_at_quad[j] = func (x_deref, y_deref);
    }

    // Project function onto DG basis using Dunavant quadrature:
    // ∫ f(x) φᵢ(x) dx ≈ Σ w_j f(x_j) φᵢ(x_j) * scaling
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto j = 0u; j < num_q; ++j)
        for (auto k = 0u; k < Base::U_DIM; ++k)
          sum[k] += Base::basis.quad.weights[j] * f_at_quad[j][k] * scaling_factor * basis_at_quad[j][i];

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
  // Reconstruction (inverse of projection)
  //=============================================================================

  /**
   * @brief Evaluate the reconstructed DG solution at a physical point
   *
   * Inverse of project_impl: reorders the vertices (data.order), maps the
   * physical point to the reference triangle and sums the orthonormal Dubiner
   * modes scaled by normalization(vol).
   *
   * @param tree_idx Tree index in forest
   * @param element Pointer to element (supplies the cell geometry)
   * @param data Leaf data holding the DG coefficients, vertex order and volume
   * @param x_phys Physical evaluation point
   * @return Solution value per component
   */
  std::array<double, Base::U_DIM>
  evaluate (int tree_idx, const t8_element_t *element, const element_t &data,
            const std::array<double, Base::DIM> &x_phys)
  {
    double vertices[3][3];
    element_vertex_coords (tree_idx, element, data.order, vertices);

    auto [trafo_mat, perm] = Base::basis.trafo_matrix_to_ref_element (vertices);
    const auto ref = Base::basis.ref_point (trafo_mat, perm, { x_phys[0], x_phys[1] });

    // ref is barycentric {lambda0, lambda1, lambda2}; the reference triangle
    // coordinates are (tau1, tau2) = (lambda1, lambda2).
    return Base::evaluate_reference (data, { ref[1], ref[2] });
  }

  /**
   * @brief Evaluate the solution gradient at a physical point
   *
   * grad[u][d] = d(u_u)/d(x_d). The basis lives in barycentric coordinates
   * (lambda0, lambda1); the chain rule uses dlambda/dx from the same transform
   * that maps a physical point to reference.
   */
  std::array<std::array<double, Base::DIM>, Base::U_DIM>
  evaluate_gradient (int tree_idx, const t8_element_t *element, const element_t &data,
                     const std::array<double, Base::DIM> &x_phys)
  {
    double vertices[3][3];
    element_vertex_coords (tree_idx, element, data.order, vertices);

    auto [trafo_mat, perm] = Base::basis.trafo_matrix_to_ref_element (vertices);
    const auto ref = Base::basis.ref_point (trafo_mat, perm, { x_phys[0], x_phys[1] });
    const auto ref_grad = Base::basis.basis_gradient ({ ref[0], ref[1] });
    const auto scaling = basis<element_t::Shape, Base::P_DIM>::normalization (data.vol);

    // dlambda/dx and dlambda/dy solve trafo * dlambda = e_x / e_y.
    std::vector<double> dlambda_dx = { 1.0, 0.0, 0.0 };
    std::vector<double> dlambda_dy = { 0.0, 1.0, 0.0 };
    lu_solve (trafo_mat, perm, dlambda_dx);
    lu_solve (trafo_mat, perm, dlambda_dy);

    std::array<std::array<double, Base::DIM>, Base::U_DIM> grad = {};
    for (auto u = 0u; u < Base::U_DIM; ++u)
      for (auto i = 0u; i < Base::DOF; ++i) {
        const double c = data.u_coeffs[element_t::dg_idx (u, i)] * scaling;
        const double dphi_dl0 = ref_grad[0][i];
        const double dphi_dl1 = ref_grad[1][i];
        grad[u][0] += c * (dphi_dl0 * dlambda_dx[0] + dphi_dl1 * dlambda_dx[1]);
        grad[u][1] += c * (dphi_dl0 * dlambda_dy[0] + dphi_dl1 * dlambda_dy[1]);
      }

    return grad;
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
  post_adaptation_hook ()
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
