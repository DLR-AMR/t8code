#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/base.hxx"
#include "t8_mra/core/adaptation.hxx"
#include "t8_mra/core/cell_geometry.hxx"
#include "t8_mra/num/geometry.hxx"
#include "t8_mra/num/mask_coefficients.hxx"
#include "t8_mra/num/nodal_to_modal.hxx"

#include <array>
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
  public multiscale_base<multiscale<TShape, U, P>, TShape, U, P>,
  public multiscale_adaptation<multiscale<TShape, U, P>> {

 public:
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
   * @brief Constructor for cartesian multiscale
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
   * @brief Detail 2-norm per component. Cartesian modes are reference-orthonormal
   *        (M = vol*I), so the raw coefficient 2-norm is already the volume-scaled
   *        detail measure the threshold expects (matches the triangle's).
   */
  std::array<double, Base::U_DIM>
  local_detail_norm (const levelmultiindex &lmi)
  {
    std::array<double, Base::U_DIM> detail_norm = {};
    const auto &details = Base::d_map.get (lmi).d_coeffs;

    for (auto u = 0u; u < Base::U_DIM; ++u) {
      double norm_sq = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
        for (auto i = 0u; i < Base::DOF; ++i) {
          const auto d = details[Base::detail_t::wavelet_idx (k, u, i)];
          norm_sq += d * d;
        }
      detail_norm[u] = std::sqrt (norm_sq);
    }

    return detail_norm;
  }

  //=============================================================================
  // Element-Specific Implementation: Projection
  //=============================================================================

  /// Corner coordinates of the element. QUAD corners are permuted (t8code swaps
  /// 2 and 3) so index 0 is the lower and the last the upper corner, as
  /// extract_cartesian_vertices expects.
  void
  element_vertex_coords (int tree_idx, const t8_element_t *element, double vertices[8][3])
  {
    constexpr int num_vertices = (Base::DIM == 1 ? 2 : (Base::DIM == 2 ? 4 : 8));
    for (int i = 0; i < num_vertices; ++i) {
      int v = i;
      if constexpr (Base::DIM == 2 && TShape == T8_ECLASS_QUAD) {
        constexpr int vertex_perm[4] = { 0, 1, 3, 2 };
        v = vertex_perm[i];
      }
      t8_forest_element_coordinate (Base::get_forest (), tree_idx, element, v, vertices[i]);
    }
  }

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
    double vertices[8][3] = {};
    element_vertex_coords (tree_idx, element, vertices);
    std::array<double, Base::DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<Base::DIM> (vertices, vertices_min, vertices_max);
    const auto geom = cell_geometry<TShape, Base::P_DIM>::from_box (
      vertices_min, vertices_max, t8_forest_element_volume (Base::get_forest (), tree_idx, element));

    // Orthonormal Legendre basis on the reference cell [0,1]^DIM: no volume/Jacobian
    // scaling in the projection. Reference quad points are mapped to physical by geom.
    const auto num_q = Base::basis.quad.num_points;
    std::vector<std::array<double, Base::DOF>> basis_at_quad (num_q);
    std::vector<std::array<double, Base::DIM>> phys_at_quad (num_q);
    std::array<double, Base::DIM> x_ref;
    for (auto q = 0u; q < num_q; ++q) {
      for (unsigned int d = 0; d < Base::DIM; ++d)
        x_ref[d] = Base::basis.quad.points[Base::DIM * q + d];
      basis_at_quad[q] = Base::basis.basis_value (x_ref);
      phys_at_quad[q] = geom.to_physical (x_ref);
    }

    // Project onto each basis function
    for (auto i = 0u; i < Base::DOF; ++i) {
      std::array<double, Base::U_DIM> sum = {};

      for (auto q = 0u; q < num_q; ++q) {
        const auto &x_phys = phys_at_quad[q];

        // Two supported func conventions:
        //   func(x, y, z) -> std::array<double, U>
        //   func(x, y, z, double* out)
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

        for (auto u = 0u; u < Base::U_DIM; ++u)
          sum[u] += Base::basis.quad.weights[q] * f_val[u] * basis_at_quad[q][i];
      }

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

    data.vol = t8_forest_element_volume (Base::get_forest (), tree_idx, element);
    project_impl (data.u_coeffs, tree_idx, element, func);

    return data;
  }

  //=============================================================================
  // Reconstruction (inverse of projection)
  //=============================================================================

  /**
   * @brief Evaluate the reconstructed DG solution at a physical point
   *
   * Inverse of project_impl: maps the physical point into the reference cell
   * (axis-aligned box) and sums the orthonormal Legendre modes. No volume
   * scaling (cartesian normalization is 1).
   *
   * @param tree_idx Tree index in forest
   * @param element Pointer to element (supplies the cell geometry)
   * @param data Leaf data holding the DG coefficients
   * @param x_phys Physical evaluation point
   * @return Solution value per component
   */
  std::array<double, Base::U_DIM>
  evaluate (int tree_idx, const t8_element_t *element, const element_t &data,
            const std::array<double, Base::DIM> &x_phys)
  {
    double vertices[8][3] = {};
    element_vertex_coords (tree_idx, element, vertices);

    std::array<double, Base::DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<Base::DIM> (vertices, vertices_min, vertices_max);

    const auto geom = cell_geometry<TShape, Base::P_DIM>::from_box (vertices_min, vertices_max, data.vol);
    return Base::evaluate_reference (data, geom.to_reference (x_phys));
  }

  /**
   * @brief Evaluate the solution gradient at a physical point
   *
   * grad[u][d] = d(u_u)/d(x_d). The box map has a diagonal Jacobian, so the
   * reference gradient is scaled by 1/(x_max - x_min) per axis.
   */
  std::array<std::array<double, Base::DIM>, Base::U_DIM>
  evaluate_gradient (int tree_idx, const t8_element_t *element, const element_t &data,
                     const std::array<double, Base::DIM> &x_phys)
  {
    double vertices[8][3] = {};
    element_vertex_coords (tree_idx, element, vertices);

    std::array<double, Base::DIM> vertices_min, vertices_max;
    extract_cartesian_vertices<Base::DIM> (vertices, vertices_min, vertices_max);

    const auto geom = cell_geometry<TShape, Base::P_DIM>::from_box (vertices_min, vertices_max, data.vol);
    const auto x_ref = geom.to_reference (x_phys);

    std::array<std::array<double, Base::DIM>, Base::U_DIM> grad = {};
    for (auto u = 0u; u < Base::U_DIM; ++u)
      grad[u] = geom.gradient (std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], Base::DOF), x_ref);

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
    Base::grid.forest = t8_forest_new_uniform (mesh, scheme, level, 0, Base::get_comm ());
    Base::build_lmi_map (
      [&] (int tree_idx, const t8_element_t *element) { return project_leaf (tree_idx, element, func); });
  }

  /**
   * @brief Load per-cell nodal DG values onto an existing forest as modal coeffs.
   *
   * Refs the committed forest (cleanup releases it); cell_nodal_values returns a
   * cell's U*DOF nodal values (component-major u*DOF + j) at reference `nodes`.
   */
  template <typename CellNodalValues>
  void
  initialize_data_nodal (t8_forest_t forest, const std::array<std::array<double, Base::DIM>, Base::DOF> &nodes,
                         CellNodalValues &&cell_nodal_values)
  {
    const nodal_to_modal<TShape, U, P> to_modal (nodes);

    t8_forest_ref (forest);
    Base::grid.forest = forest;

    Base::build_lmi_map ([&] (int tree_idx, const t8_element_t *element) {
      element_t data;
      data.vol = t8_forest_element_volume (Base::get_forest (), tree_idx, element);
      const auto nodal = cell_nodal_values (tree_idx, element);
      to_modal (std::span<const double> (nodal.data (), nodal.size ()), data.u_coeffs);

      return data;
    });
  }

  /**
   * @brief Reconstruct per-cell nodal DG values from the current forest.
   *
   * Inverse of initialize_data_nodal: walks the (adapted) forest in SFC order
   * and hands each leaf's nodal values at reference `nodes` to
   * write_cell_nodal_values (tree_idx, element, span<const double> of U*DOF,
   * component-major u*DOF + j).
   */
  template <typename WriteCellNodalValues>
  void
  export_data_nodal (const std::array<std::array<double, Base::DIM>, Base::DOF> &nodes,
                     WriteCellNodalValues &&write_cell_nodal_values)
  {
    const modal_to_nodal<TShape, U, P> to_nodal (nodes);
    const auto *scheme = t8_forest_get_scheme (Base::get_forest ());
    auto *lmi_map = Base::get_lmi_map ();

    Base::for_each_local_leaf (
      [&] (t8_locidx_t tree_idx, const t8_element_t *element, unsigned int, t8_gloidx_t global_tree) {
        const levelmultiindex lmi (global_tree, element, scheme);
        const auto *data = lmi_map->find (lmi);
        const auto nodal = to_nodal (std::span<const double> (data->u_coeffs.data (), data->u_coeffs.size ()));
        write_cell_nodal_values (tree_idx, element, std::span<const double> (nodal.data (), nodal.size ()));
      });
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
