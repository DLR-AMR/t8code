#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hpp"
#include "t8_mra/t8_basis.hpp"

namespace t8_mra
{

/// Specialization of multiscale_data for cartesian elements
template <t8_eclass TShape>
  requires is_cartesian<TShape>
struct multiscale_data<TShape>
{
  static constexpr unsigned short DIM = TShape == T8_ECLASS_LINE ? 1 : (TShape == T8_ECLASS_QUAD ? 2 : 3);

  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

/// Cartesian-specific multiscale class
template <t8_eclass TShape, int U, int P>
  requires is_cartesian<TShape>
class multiscale<TShape, U, P>: public multiscale_data<TShape> {
 public:
  using element_t = data_per_element<TShape, U, P>;
  using levelmultiindex = t8_mra::levelmultiindex<TShape>;
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;
  static constexpr auto Shape = TShape;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;

  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

  using multiscale_data<TShape>::mask_coefficients;
  using multiscale_data<TShape>::inverse_mask_coefficients;

 public:
  unsigned int maximum_level;
  double c_thresh;
  std::array<double, U> c_scaling;
  int gamma;
  std::vector<double> eps;
  int num_quad_points_1d;  // Number of 1D quadrature points
  t8_mra::dg_basis<element_t> DG_basis;

  levelindex_map<levelmultiindex, element_t> d_map;
  levelindex_set<levelmultiindex> td_set;
  levelindex_set<levelmultiindex> refinement_set;
  levelindex_set<levelmultiindex> coarsening_set;

  /// Forest data
  t8_forest_t forest;
  bool balanced;
  sc_MPI_Comm comm;

 public:
  /**
   * @brief Constructor for cartesian multiscale
   *
   * @param _max_level Maximum refinement level
   * @param _c_thresh Coarsening threshold
   * @param _gamma Gamma parameter for thresholding
   * @param _num_quad_points_1d Number of 1D quadrature points (Gauss-Legendre)
   * @param _balanced Whether to balance the forest
   * @param _comm MPI communicator
   */
  multiscale (int _max_level, double _c_thresh, int _gamma, int _num_quad_points_1d, bool _balanced, sc_MPI_Comm _comm)
    : maximum_level (_max_level), c_thresh (_c_thresh), gamma (_gamma), num_quad_points_1d (_num_quad_points_1d),
      balanced (_balanced), comm (_comm), DG_basis (_num_quad_points_1d, P_DIM), d_map (maximum_level),
      td_set (maximum_level), refinement_set (maximum_level), coarsening_set (maximum_level)
  {
    // TODO: Initialize mask coefficients for cartesian elements
    // For now, we skip this - it's needed for multiscale transformation
    // t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, mask_coefficients, inverse_mask_coefficients);
  }

  t8_forest_t
  get_forest ()
  {
    return forest;
  }

  t8_mra::forest_data<element_t> *
  get_user_data ()
  {
    return reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest));
  }

  t8_mra::levelindex_map<levelmultiindex, element_t> *
  get_lmi_map ()
  {
    return get_user_data ()->lmi_map;
  }

  /**
   * @brief Projects a function onto the DG basis for a cartesian element
   *
   * For cartesian elements (LINE, QUAD, HEX), vertices define an axis-aligned box.
   * The projection uses Gauss-Legendre quadrature on the reference element [0,1]^DIM.
   *
   * @tparam Func Function type with signature std::array<double, U_DIM>(double, double, ...)
   * @param dg_coeffs Output vector for DG coefficients
   * @param tree_idx Tree index in forest
   * @param element Pointer to element
   * @param func Function to project (takes DIM coordinates, returns U_DIM values)
   */
  template <typename Func>
  void
  project (std::vector<double> &dg_coeffs, int tree_idx, const t8_element_t *element, Func &&func)
  {
    // Extract element vertices
    constexpr int num_vertices = (DIM == 1 ? 2 : (DIM == 2 ? 4 : 8));
    double vertices[8][3] = {};  // Max 8 vertices for HEX

    for (int i = 0; i < num_vertices; ++i) {
      t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[i]);
    }

    // Get physical quadrature points
    const auto phys_quad_points = DG_basis.deref_quad_points (vertices);

    // DEBUG: Print first few physical quad points for first element
    static bool printed_debug = false;
    if (!printed_debug && tree_idx == 0) {
      printf ("DEBUG: Physical quadrature points for first element:\n");
      for (int i = 0; i < std::min ((size_t) 4u, DG_basis.num_quad_points); ++i) {
        printf ("  quad[%d]: (%.6f, %.6f)\n", i, phys_quad_points[2 * i], phys_quad_points[2 * i + 1]);
      }
      printf ("DEBUG: Reference quadrature points:\n");
      for (int i = 0; i < std::min ((size_t) 4u, DG_basis.num_quad_points); ++i) {
        printf ("  ref[%d]: (%.6f, %.6f)\n", i, DG_basis.ref_quad_points[2 * i], DG_basis.ref_quad_points[2 * i + 1]);
      }
      printed_debug = true;
    }

    // Get element volume (Jacobian determinant for coordinate transform)
    const auto volume = t8_forest_element_volume (forest, tree_idx, element);

    // Evaluate all basis functions at all reference quadrature points (precompute)
    std::vector<std::array<double, DOF>> basis_at_quad (DG_basis.num_quad_points);
    for (auto j = 0u; j < DG_basis.num_quad_points; ++j) {
      // Extract reference quadrature point
      std::vector<double> x_ref (DIM);
      for (unsigned int d = 0; d < DIM; ++d) {
        x_ref[d] = DG_basis.ref_quad_points[DIM * j + d];
      }
      const auto basis_vals = DG_basis.basis_value (x_ref);
      for (auto i = 0u; i < DOF; ++i) {
        basis_at_quad[j][i] = basis_vals[i];
      }
    }

    // Project onto each basis function
    for (auto i = 0u; i < DOF; ++i) {
      std::array<double, U_DIM> sum = {};

      for (auto j = 0u; j < DG_basis.num_quad_points; ++j) {
        // Extract physical coordinates from flattened array
        std::array<double, DIM> x_phys;
        for (unsigned int d = 0; d < DIM; ++d) {
          x_phys[d] = phys_quad_points[DIM * j + d];
        }

        // Evaluate function at physical quadrature point
        std::array<double, U_DIM> f_val;
        if constexpr (DIM == 1) {
          f_val = func (x_phys[0]);
        }
        else if constexpr (DIM == 2) {
          f_val = func (x_phys[0], x_phys[1]);
        }
        else if constexpr (DIM == 3) {
          f_val = func (x_phys[0], x_phys[1], x_phys[2]);
        }

        // Accumulate quadrature sum: integral(f * phi_i)
        // Note: For orthonormal basis, the volume scaling cancels out
        for (auto k = 0u; k < U_DIM; ++k) {
          sum[k] += DG_basis.quad_weights[j] * f_val[k] * basis_at_quad[j][i];
        }
      }

      // Store coefficients directly (already includes volume scaling from quadrature)
      for (auto k = 0u; k < U_DIM; ++k) {
        dg_coeffs[element_t::dg_idx (k, i)] = sum[k];
      }
    }
  }

  /**
   * @brief Initializes uniform forest with projected data
   *
   * Creates a uniform forest at the specified level and projects the given function
   * onto all elements.
   *
   * @tparam Func Function type
   * @param mesh Coarse mesh
   * @param scheme Element scheme
   * @param level Uniform refinement level
   * @param func Function to project
   */
  template <typename Func>
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, Func &&func)
  {
    forest = t8_forest_new_uniform (mesh, scheme, level, 0, comm);

    levelmultiindex *elem_data;
    t8_mra::forest_data<element_t> *user_data;

    user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    elem_data = T8_ALLOC (levelmultiindex, 1);

    T8_ASSERT (t8_forest_is_commited (forest));

    const auto num_local_elements = t8_forest_get_global_num_leaf_elements (forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (forest);

    user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;

    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto base_element = t8_forest_global_tree_id (forest, tree_idx);

      for (auto ele_idx = 0u; ele_idx < num_elements_in_tree; ++ele_idx, ++current_idx) {
        element_t data_element;
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (base_element, element, scheme);

        // Cartesian elements don't need vertex ordering
        data_element.vol = t8_forest_element_volume (forest, tree_idx, element);

        // Project function onto element
        project (data_element.u_coeffs, tree_idx, element, func);

        user_data->lmi_map->insert (lmi, data_element);

        // Insert lmi into forest
        t8_mra::set_lmi_forest_data (user_data, current_idx, lmi);
      }
    }

    T8_FREE (elem_data);
    t8_forest_set_user_data (forest, user_data);
  }

  /**
   * @brief Evaluates the solution at a physical point
   *
   * @param tree_idx Tree index
   * @param element Element pointer
   * @param x Physical coordinates (DIM values)
   * @return std::array<double, U_DIM> Solution values at point x
   */
  std::array<double, U_DIM>
  evaluate_at_point (int tree_idx, const t8_element_t *element, const std::vector<double> &x)
  {
    // Extract element vertices
    constexpr int num_vertices = (DIM == 1 ? 2 : (DIM == 2 ? 4 : 8));
    double vertices[8][3] = {};

    for (int i = 0; i < num_vertices; ++i) {
      t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[i]);
    }

    // Map physical point to reference element
    const auto x_ref = DG_basis.ref_point (vertices, x);

    // Evaluate basis functions at reference point
    const auto basis_val = DG_basis.basis_value (x_ref);

    // Get element's LMI and coefficients
    const auto offset = t8_forest_get_tree_element_offset (forest, tree_idx);
    const auto elem_idx = offset;  // TODO: Need actual element index
    const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);
    const auto &coeffs = get_user_data ()->lmi_map->get (lmi).u_coeffs;

    // Reconstruct solution
    std::array<double, U_DIM> result = {};
    for (auto k = 0u; k < U_DIM; ++k) {
      for (auto i = 0u; i < DOF; ++i) {
        result[k] += coeffs[element_t::dg_idx (k, i)] * basis_val[i];
      }
    }

    return result;
  }

  void
  cleanup ()
  {
    delete get_user_data ()->lmi_map;
    sc_array_destroy (get_user_data ()->lmi_idx);
    T8_FREE (get_user_data ());
    t8_forest_unref (&forest);
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
