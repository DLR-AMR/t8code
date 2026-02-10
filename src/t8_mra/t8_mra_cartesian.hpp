#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra.hpp"
#include "t8_mra/t8_basis.hpp"
#include "t8_mra/num/mask_coefficients.hpp"
#include "t8_mra/num/mask_coefficients_compute.hxx"

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
    // Use COMPUTED mask coefficients to ensure correctness for linear functions
    std::cout << "=== Initializing Mask Coefficients (COMPUTED) ===\n";
    t8_mra::initialize_mask_coefficients_computed<TShape> (P_DIM, DOF, mask_coefficients, inverse_mask_coefficients);
    // t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, mask_coefficients, inverse_mask_coefficients);

    // Debug: Print mask coefficient M0[0] (first child) to verify prolongation property
    std::cout << "\n=== DEBUG: Mask Coefficients M0[0] (child 0) ===\n";
    std::cout << "For constant f=1, M0[k]*[1,0,...,0] should give [1,0,...,0]\n";
    for (auto k = 0u; k < 4; ++k) {
      std::cout << "For linear, M" << k << "[0] should account for shift\n\n";
      std::cout << "M" << k << "[0] matrix (" << DOF << "x" << DOF << "):\n";
      for (size_t i = 0; i < DOF; ++i) {
        std::cout << "  Row " << i << ": ";
        for (size_t j = 0; j < DOF; ++j) {
          std::cout << std::setw (12) << std::setprecision (6) << std::fixed << mask_coefficients[k](i, j) << " ";
        }
        std::cout << "\n";
      }
    }

    // Test prolongation of constant function: M0[0][:,0] should be [1, ?, ?, ...]
    std::cout << "\nTest prolongation of constant [1,0,0,...]: M0[0] * [1, 0, 0, ...] = \n  [";
    for (size_t i = 0; i < DOF; ++i) {
      std::cout << std::setw (12) << std::setprecision (6) << std::fixed << mask_coefficients[0](i, 0);
      if (i < DOF - 1)
        std::cout << ", ";
    }
    std::cout << "]\nExpected: [1, 0, 0, ...] for constant to be preserved exactly\n";

    // Check all 4 children's first column (constant mode prolongation)
    std::cout << "\nAll children M[k][:,0] (constant mode prolongation):\n";
    for (size_t k = 0; k < levelmultiindex::NUM_CHILDREN; ++k) {
      std::cout << "  Child " << k << ": [";
      const size_t max_print = (DOF < 4) ? DOF : 4;
      for (size_t i = 0; i < max_print; ++i) {
        std::cout << std::setw (10) << std::setprecision (4) << std::fixed << mask_coefficients[k](i, 0);
        if (i < max_print - 1)
          std::cout << ", ";
      }
      if (DOF > 4)
        std::cout << ", ...";
      std::cout << "]\n";
    }
    std::cout << "\n";
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

    if constexpr (DIM == 2 && TShape == T8_ECLASS_QUAD) {
      // t8code QUAD vertex ordering: 0-1-2-3 as (0,0)-(1,0)-(0,1)-(1,1)
      // But we need: 0-1-2-3 as (0,0)-(1,0)-(1,1)-(0,1) for standard quad
      // So we need to swap vertices 2 and 3
      const int vertex_perm[4] = { 0, 1, 3, 2 };
      for (int i = 0; i < num_vertices; ++i) {
        t8_forest_element_coordinate (forest, tree_idx, element, vertex_perm[i], vertices[i]);
        // std::cout << "[" << vertices[i][0] << " " << vertices[i][1] << "]";
      }
      // std::cout << "\n";
    }
    else {
      for (int i = 0; i < num_vertices; ++i) {
        t8_forest_element_coordinate (forest, tree_idx, element, i, vertices[i]);
      }
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
        for (auto u = 0u; u < U_DIM; ++u) {
          sum[u] += DG_basis.quad_weights[j] * f_val[u] * basis_at_quad[j][i];
        }
      }

      // Store coefficients directly (already includes volume scaling from quadrature)
      for (auto u = 0u; u < U_DIM; ++u) {
        dg_coeffs[element_t::dg_idx (u, i)] = sum[u];
        std::cout << std::setprecision (6) << sum[u] << " ";
      }
    }
    std::cout << std::endl;
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

    T8_ASSERT (t8_forest_is_committed (forest));

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

    // Debug: Print first few LMIs to verify they're stored correctly
    std::cout << "=== DEBUG: First 5 LMIs in map ===\n";
    int count = 0;
    for (const auto &[lmi, data] : user_data->lmi_map->operator[] (level)) {
      std::cout << "  LMI " << count << ": index=" << lmi.index << ", level=" << lmi.level () << "\n";
      if (++count >= 5)
        break;
    }
    std::cout << "Total elements in lmi_map at level " << level << ": "
              << user_data->lmi_map->operator[] (level).size () << "\n";
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

  /**
   * @brief Multiscale transformation: compute parent coefficients and detail coefficients from children
   *
   * For cartesian elements, this is simpler than triangles as there's no vertex reordering needed.
   *
   * @param l_min Minimum level
   * @param l_max Maximum level
   */
  // void
  // multiscale_transformation (unsigned int l_min, unsigned int l_max)
  // {
  //   index_set I_set;
  //   element_t data_on_coarse;
  //   std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;
  //   bool first_element_printed = false;
  //
  //   // Match triangle implementation: no explicit mask_factor
  //   // Normalization is baked into the mask coefficients
  //   // const double mask_factor = 1.0;
  //   const double mask_factor = 0.5;
  //   const double inv_mask_factor = 1.0;
  //
  //   for (auto l = l_max; l > l_min; --l) {
  //     for (const auto &[lmi, _] : get_user_data ()->lmi_map->operator[] (l))
  //       I_set.emplace (t8_mra::parent_lmi (lmi));
  //
  //     d_map[l - 1].reserve (get_user_data ()->lmi_map->size (l));
  //
  //     for (const auto &lmi : I_set) {
  //       const auto siblings_lmi = t8_mra::children_lmi (lmi);
  //
  //       // Load children in mask coefficient order
  //       for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
  //         data_on_siblings[k] = get_user_data ()->lmi_map->get (siblings_lmi[k]);
  //       }
  //
  //       for (auto u = 0u; u < U_DIM; ++u) {
  //         /// Single scale of lmi: compute parent coefficients from children with normalization
  //         for (auto i = 0u; i < DOF; ++i) {
  //           auto sum = 0.0;
  //
  //           for (auto j = 0u; j < DOF; ++j)
  //             for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
  //               // Match triangle implementation: use (j,i) for forward MST
  //               // sum
  //               // += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i) * mask_factor;
  //               sum
  //                 += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j) * mask_factor;
  //
  //           // Debug: print forward MST calculation for coeff 1
  //           static bool fwd_mst_printed = false;
  //           if (!fwd_mst_printed && l == l_max && u == 0 && i == 1) {
  //             fwd_mst_printed = true;
  //             std::cout << "\n=== DEBUG: Forward MST for u_parent[1] ===\n";
  //             std::cout << "Child coefficients u_child[1]: [";
  //             for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
  //               std::cout << data_on_siblings[k].u_coeffs[element_t::dg_idx (u, 1)]
  //                         << (k < levelmultiindex::NUM_CHILDREN - 1 ? ", " : "");
  //             std::cout << "]\n";
  //             std::cout << "Mask matrices M_k[:,1] (column 1 of each mask):\n";
  //             for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
  //               std::cout << "  M_" << k << "[:,1] = [";
  //               for (auto j = 0u; j < DOF; ++j)
  //                 std::cout << mask_coefficients[k](j, 1) << (j < DOF - 1 ? ", " : "");
  //               std::cout << "]\n";
  //             }
  //             std::cout << "Computing: u_parent[1] = sum_k sum_j M_k[j,1] * u_child_k[j]\n";
  //             for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
  //               double contrib = 0.0;
  //               for (auto j = 0u; j < DOF; ++j)
  //                 contrib += mask_coefficients[k](j, 1) * data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)];
  //               std::cout << "  Child " << k << " contribution: " << contrib << "\n";
  //             }
  //             std::cout << "Sum before mask_factor: " << sum << "\n";
  //             std::cout << "u_parent[1] = " << (sum * mask_factor) << "\n\n";
  //           }
  //
  //           // data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = sum * mask_factor;  //! Check that
  //           data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = sum;  //! Check that
  //           data_on_coarse.vol = data_on_siblings[0].vol * levelmultiindex::NUM_CHILDREN;
  //         }
  //
  //         /// Details as differences: d_k = u_k - M_k * u_parent
  //         /// Note: u_parent was scaled by mask_factor, so we need to divide it back out
  //         for (auto i = 0u; i < DOF; ++i)
  //           for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
  //             auto sum = 0.0;
  //             for (auto j = 0u; j < DOF; ++j)
  //               // sum += mask_coefficients[k](i, j) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];
  //               sum += mask_coefficients[k](j, i) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];
  //
  //             data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, i)]
  //               = data_on_siblings[k].u_coeffs[element_t::dg_idx (u, i)] - sum;
  //
  //             // Debug: print detail calculation for coeff 1
  //             const auto detail = data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, i)];
  //             static bool detail_calc_printed = false;
  //             if (!detail_calc_printed && l == l_max && u == 0 && k == 0 && i == 1) {
  //               detail_calc_printed = true;
  //               std::cout << "\n=== DEBUG: Detail calculation for child 0, coeff 1 ===\n";
  //               std::cout << "u_parent coeffs: [";
  //               for (auto jj = 0u; jj < DOF; ++jj)
  //                 std::cout << data_on_coarse.u_coeffs[element_t::dg_idx (u, jj)] << (jj < DOF - 1 ? ", " : "");
  //               std::cout << "]\n";
  //               std::cout << "M_0[1,:] (row 1 of mask for child 0): [";
  //               for (auto jj = 0u; jj < DOF; ++jj)
  //                 std::cout << mask_coefficients[0](1, jj) << (jj < DOF - 1 ? ", " : "");
  //               std::cout << "]\n";
  //               std::cout << "Computing M_0[1,:] * u_parent:\n";
  //               double check_sum = 0.0;
  //               for (auto jj = 0u; jj < DOF; ++jj) {
  //                 double contrib = mask_coefficients[0](1, jj) * data_on_coarse.u_coeffs[element_t::dg_idx (u, jj)];
  //                 std::cout << "  j=" << jj << ": M[1," << jj << "]=" << mask_coefficients[0](1, jj) << " * u_parent["
  //                           << jj << "]=" << data_on_coarse.u_coeffs[element_t::dg_idx (u, jj)] << " = " << contrib
  //                           << "\n";
  //                 check_sum += contrib;
  //               }
  //               std::cout << "Sum (M*u_parent): " << check_sum << "\n";
  //               std::cout << "u_child[1] = " << data_on_siblings[k].u_coeffs[element_t::dg_idx (u, i)] << "\n";
  //               std::cout << "detail = " << detail << "\n\n";
  //             }
  //
  //             // Debug: check if children are identical (constant function test)
  //             static bool detail_debug = false;
  //             if (!detail_debug && u == 0 && i == 0 && k < 4) {
  //               bool all_same = true;
  //               for (auto kk = 1u; kk < 4; ++kk) {
  //                 if (std::abs (data_on_siblings[kk].u_coeffs[0] - data_on_siblings[0].u_coeffs[0]) > 1e-10) {
  //                   all_same = false;
  //                   break;
  //                 }
  //               }
  //               if (all_same) {
  //                 std::cout << "        DEBUG: All 4 children have same u_coeffs[0] = "
  //                           << data_on_siblings[0].u_coeffs[0] << "\n";
  //                 std::cout << "          Parent u_coeffs[0] = " << data_on_coarse.u_coeffs[0] << "\n";
  //                 std::cout << "          M_" << k << " * u_parent[0] = " << sum << "\n";
  //                 std::cout << "          Detail d_" << k << "[0] = " << (data_on_siblings[k].u_coeffs[0] - sum)
  //                           << "\n";
  //                 if (k == 3)
  //                   detail_debug = true;
  //               }
  //             }
  //           }
  //       }
  //
  //       // Cartesian elements don't need vertex order updates
  //       get_user_data ()->lmi_map->insert (lmi, data_on_coarse);
  //       d_map[l - 1].emplace (lmi, data_on_coarse);
  //     }
  //
  //     get_user_data ()->lmi_map->erase (l);
  //     I_set.clear ();
  //   }
  // }
  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    index_set I_set;
    element_t data_on_coarse;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;

    const double mask_factor = 1.0;
    const double inv_mask_factor = 1.0;

    for (auto l = l_max; l > l_min; --l) {
      for (const auto &[lmi, _] : get_user_data ()->lmi_map->operator[] (l))
        I_set.emplace (t8_mra::parent_lmi (lmi));

      d_map[l - 1].reserve (get_user_data ()->lmi_map->size (l));

      for (const auto &lmi : I_set) {
        const auto siblings_lmi = t8_mra::children_lmi (lmi);

        // Load children in mask coefficient order
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          data_on_siblings[k] = get_user_data ()->lmi_map->get (siblings_lmi[k]);

        // Coarser u_coeff
        for (auto u = 0u; u < U_DIM; ++u) {
          for (auto i = 0u; i < DOF; ++i) {
            auto sum = 0.0;
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
              for (auto j = 0u; j < DOF; ++j) {
                sum += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i);
                // sum += data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j);
              }
            }
            data_on_coarse.vol = data_on_siblings[0].vol * levelmultiindex::NUM_CHILDREN;
            data_on_coarse.u_coeffs[element_t::dg_idx (u, i)] = sum / levelmultiindex::NUM_CHILDREN;

            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
              for (auto j = 0u; j < DOF; ++j)
                data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, j)]
                  = data_on_siblings[k].u_coeffs[element_t::dg_idx (u, j)];
          }
        }

        for (auto u = 0u; u < U_DIM; ++u) {
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
            for (auto i = 0u; i < DOF; ++i) {
              auto sum = 0.0;
              for (auto j = 0u; j < DOF; ++j) {
                sum += mask_coefficients[k](i, j) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];
                // sum += mask_coefficients[k](j, i) * data_on_coarse.u_coeffs[element_t::dg_idx (u, j)];
              }
              data_on_coarse.d_coeffs[element_t::wavelet_idx (k, u, i)] -= sum;
            }
          }
        }

        for (auto i = 0u; i < DOF; ++i) {
          std::cout << std::scientific << "u: " << data_on_coarse.u_coeffs[element_t::dg_idx (0, i)] << " ";
        }
        std::cout << "\n";

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          std::cout << "d" << k << ": ";
          for (auto i = 0u; i < DOF; ++i) {
            std::cout << std::scientific << data_on_coarse.d_coeffs[element_t::wavelet_idx (k, 0, i)] << " ";
          }
          std::cout << "\n";
        }
        std::cout << "\n";

        get_user_data ()->lmi_map->insert (lmi, data_on_coarse);
        d_map[l - 1].emplace (lmi, data_on_coarse);
      }

      get_user_data ()->lmi_map->erase (l);

      I_set.clear ();
    }
  }

  /**
   * @brief Inverse multiscale transformation: reconstruct children from parent + details
   *
   * @param l_min Minimum level
   * @param l_max Maximum level
   */
  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    element_t new_data;
    std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_children;
    static bool debug_printed = false;

    // Our u_coeffs are L²-normalized physical coefficients
    // Inverse MST: u_child_k = M_k * u_parent + d_k (no additional scaling)
    const double inv_mask_factor = 1.0;
    // const double inv_mask_factor = 0.5;

    for (auto l = l_min; l < l_max; ++l) {
      get_user_data ()->lmi_map->operator[] (l + 1).reserve (d_map[l].size ());

      std::cout << "d_map[" << l << "] size: " << d_map[l].size () << "\n";

      for (const auto &[lmi, d] : d_map[l]) {
        const auto children_lmi = t8_mra::children_lmi (lmi);
        const auto lmi_data = get_user_data ()->lmi_map->get (lmi);
        const auto &details = d.d_coeffs;

        // children_lmi already returns children in mask coefficient order
        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
          for (auto u = 0u; u < U_DIM; ++u) {
            for (auto i = 0u; i < DOF; ++i) {
              auto sum = 0.0;

              for (auto j = 0u; j < DOF; ++j)
                sum += lmi_data.u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](i, j);
              // sum += lmi_data.u_coeffs[element_t::dg_idx (u, j)] * mask_coefficients[k](j, i);

              new_data.u_coeffs[element_t::dg_idx (u, i)]
                = details[element_t::wavelet_idx (k, u, i)] + sum * inv_mask_factor;
            }
          }

          new_data.vol = lmi_data.vol / levelmultiindex::NUM_CHILDREN;
          get_user_data ()->lmi_map->insert (children_lmi[k], new_data);

          // std::cout << "DG after:\n";
          // for (auto i = 0u; i < DOF; ++i)
          //   std::cout << new_data.u_coeffs[element_t::dg_idx (0, i)] << " ";
          // std::cout << "\n";
        }
        get_user_data ()->lmi_map->erase (lmi);
      }
      d_map.erase (l);
    }
  }

  /**
   * @brief Compute local detail norm for hard thresholding
   *
   * @param lmi Level multi-index
   * @return std::array<double, U_DIM> Norm for each component
   */
  std::array<double, U_DIM>
  local_detail_norm (const levelmultiindex &lmi)
  {
    std::array<double, U_DIM> detail_norm = {};
    const auto vol = d_map.get (lmi).vol;

    // if (!d_map.contains (lmi))
    //   return detail_norm;

    const auto &data = d_map.get (lmi);

    for (auto u = 0u; u < U_DIM; ++u) {
      double sum = 0.0;
      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        for (auto i = 0u; i < DOF; ++i) {
          const auto d_val = data.d_coeffs[element_t::wavelet_idx (k, u, i)];
          sum += d_val * d_val;
        }
      }
      detail_norm[u] = std::sqrt (sum / vol);
    }

    return detail_norm;
    // std::array<double, U_DIM> tmp = {};
    // const auto vol = d_map.get (lmi).vol;
    //
    // for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
    //   for (auto u = 0u; u < U_DIM; ++u)
    //     for (auto i = 0u; i < DOF; ++i) {
    //       const auto d = d_map.get (lmi).d_coeffs[element_t::wavelet_idx (k, u, i)];
    //       tmp[u] += d * d;
    //     }
    // }
    //
    // for (auto u = 0u; u < U_DIM; ++u)
    //   tmp[u] = std::sqrt (tmp[u] / vol);
    //
    // return tmp;
  }

  /**
   * @brief Compute local threshold value based on element level and volume
   *
   * Implements uniform subdivision threshold from Veli eq. (2.44)
   *
   * @param lmi Level multi-index
   * @return double Threshold value
   */
  double
  local_threshold_value (const levelmultiindex &lmi)
  {
    const auto vol = d_map.get (lmi).vol;

    const auto level_diff = maximum_level - lmi.level ();
    const auto h_lambda = std::sqrt (vol);
    const auto h_max_level = std::pow (vol / std::pow (levelmultiindex::NUM_CHILDREN, level_diff), (gamma + 1.0) / 2.0);

    static bool debug_printed = false;
    if (!debug_printed && lmi.level () == 5) {
      std::cout << "        DEBUG threshold calc at level 5:\n";
      std::cout << "          vol = " << vol << ", level_diff = " << level_diff << "\n";
      std::cout << "          h_lambda = " << h_lambda << ", h_max_level = " << h_max_level << "\n";
      std::cout << "          threshold = " << (h_max_level / h_lambda) << "\n";
      debug_printed = true;
    }

    return h_max_level / h_lambda;
  }

  /**
   * @brief Compute threshold scaling factor
   *
   * @return std::array<double, U_DIM> Scaling factor for each component
   */
  std::array<double, U_DIM>
  threshold_scaling_factor ()
  {
    std::array<double, U_DIM> scaling;

    // Compute scaling based on (2.39) in wavelet theory
    for (auto u = 0u; u < U_DIM; ++u) {
      double max_val = 0.0;

      // Iterate over all levels
      for (auto l = 0u; l <= maximum_level; ++l) {
        for (const auto &[lmi, data] : d_map[l]) {
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
            for (auto i = 0u; i < DOF; ++i) {
              max_val = std::max (max_val, std::abs (data.d_coeffs[element_t::wavelet_idx (k, u, i)]));
            }
          }
        }
      }

      scaling[u] = max_val > 0.0 ? max_val : 1.0;
    }

    return scaling;
  }

  /**
   * @brief Hard thresholding: mark elements for coarsening if details are small
   *
   * @param l_min Minimum level
   * @param l_max Maximum level
   */
  void
  hard_thresholding (int l_min, int l_max)
  {
    int num_kept = 0;
    int num_total = 0;
    double max_norm = 0.0;
    double min_norm = 1e100;
    double max_eps = 0.0;

    for (auto l = l_min; l < l_max; ++l) {
      for (const auto &[lmi, d] : d_map[l]) {
        auto tmp = local_detail_norm (lmi);

        for (auto u = 0u; u < U_DIM; ++u)
          tmp[u] /= c_scaling[u];

        const auto local_norm = *std::max_element (tmp.begin (), tmp.end ());
        const auto local_eps = c_thresh * local_threshold_value (lmi);

        max_norm = std::max (max_norm, local_norm);
        min_norm = std::min (min_norm, local_norm);
        max_eps = std::max (max_eps, local_eps);
        num_total++;

        if (local_norm > local_eps) {
          td_set.insert (lmi);
          num_kept++;
        }
      }
    }

    if (num_total > 0) {
      std::cout << "      Thresholding: kept " << num_kept << "/" << num_total << ", detail norm range: [" << min_norm
                << ", " << max_norm << "], max threshold: " << max_eps << "\n";
    }
  }

  /**
   * @brief Generate td_tree by including parents of elements in td_set
   *
   * @param l_min Minimum level
   * @param l_max Maximum level
   */
  void
  generate_td_tree (int l_min, int l_max)
  {
    for (auto l = l_max; l > l_min; --l) {
      for (const auto &lmi : td_set[l]) {
        const auto parent_lmi = t8_mra::parent_lmi (lmi);
        if (!td_set.contains (parent_lmi))
          td_set.insert (parent_lmi);
      }
    }
  }

  /**
   * @brief Synchronize d_map with td_set: keep only elements in td_set
   *
   * @param l_min Minimum level
   * @param l_max Maximum level
   */
  void
  sync_d_with_td (int l_min, int l_max)
  {
    // Create a copy of d_map to iterate over
    auto d_copy = d_map;

    for (auto l = static_cast<int> (l_max) - 1; l >= static_cast<int> (l_min); --l) {
      for (const auto &[lmi, _] : d_copy[l]) {
        if (!td_set.contains (lmi)) {
          d_map.erase (lmi);
          for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
            coarsening_set.insert (t8_mra::children_lmi (lmi)[k]);
          }
        }
      }
    }
  }

  /**
   * @brief Coarsening callback for t8code adapt
   *
   * @return -1 to coarsen, 0 to keep
   */
  int
  coarsening_callback_new (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                           t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family,
                           const int num_elements, t8_element_t *elements[])
  {
    if (!is_family)
      return 0;

    const auto element_level = scheme->element_get_level (tree_class, elements[0]);

    if (element_level != get_user_data ()->current_refinement_level)
      return 0;

    // Get the parent LMI from the first child
    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
    const auto elem_idx = local_ele_idx + offset;
    const auto first_child_lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), elem_idx);
    const auto parent = t8_mra::parent_lmi (first_child_lmi);

    // Check if ALL children are in coarsening_set
    const auto children = t8_mra::children_lmi (parent);
    for (const auto &child : children) {
      if (!coarsening_set.contains (child))
        return 0;  // At least one child should not be coarsened
    }

    return -1;  // All children can be coarsened
  }

  /**
   * @brief Iterate-replace callback for updating user data during adaptation
   */
  void
  iterate_replace_callback_new (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                                t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    auto *new_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_new));
    auto *old_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_old));

    if (refine == -1) {
      // Coarsening: first_incoming is the parent
      // Get parent LMI from first child in old forest
      const auto first_child_lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing);
      const auto parent_lmi_from_child = t8_mra::parent_lmi (first_child_lmi);

      // Also compute from element to verify they match
      const auto *parent_elem = t8_forest_get_leaf_element_in_tree (forest_new, which_tree, first_incoming);
      const auto parent_lmi_from_elem = levelmultiindex (which_tree, parent_elem, scheme);

      // Debug: check if they match
      if (parent_lmi_from_child.index != parent_lmi_from_elem.index) {
        std::cout << "WARNING: Parent LMI mismatch!\n";
        std::cout << "  From child: " << parent_lmi_from_child.index << "\n";
        std::cout << "  From elem:  " << parent_lmi_from_elem.index << "\n";
      }

      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, parent_lmi_from_child);
    }
    else if (refine == 0) {
      // No change: copy LMI
      for (int i = 0; i < num_incoming; ++i) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing + i);
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, lmi);
      }
    }
    else if (refine == 1) {
      // Refinement: compute children LMIs
      const auto *parent_elem = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing);
      const auto parent_lmi = levelmultiindex (which_tree, parent_elem, scheme);
      const auto children = t8_mra::children_lmi (parent_lmi);

      for (int i = 0; i < num_incoming; ++i)
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, children[i]);
    }
  }

  /**
   * @brief Main coarsening function: performs multiscale decomposition and coarsening
   *
   * @param min_level Minimum refinement level
   * @param max_level Maximum refinement level
   */
  void
  coarsening_new (int min_level, int max_level)
  {
    static auto static_coarsening_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family, const int num_elements,
                t8_element_t *elements[]) -> int {
      return coarsening_callback_new (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                                      num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      iterate_replace_callback_new (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                    first_outgoing, num_incoming, first_incoming);
    };

    /// Scaling due to (2.39)
    c_scaling = threshold_scaling_factor ();

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (forest);

      get_user_data ()->current_refinement_level = l;

      // Forward MST: level -> level-1
      multiscale_transformation (l - 1, l);

      // Thresholding
      hard_thresholding (l - 1, l);
      generate_td_tree (l - 1, l);
      sync_d_with_td (min_level, max_level);

      // Inverse MST: reconstruct children for kept parents
      inverse_multiscale_transformation (l - 1, l);

      std::cout << "  Level " << l << " -> " << (l - 1) << ": kept " << td_set[l - 1].size () << "/"
                << d_map[l - 1].size () << " parents, will coarsen " << coarsening_set[l].size () << " children\n";

      new_forest = t8_forest_new_adapt (
        forest,
        [] (auto *forest, auto *forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto *scheme,
            const auto is_family, const auto num_elements, auto *elements[]) -> int {
          return static_coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, get_user_data ());

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (maximum_level);
      std::swap (new_user_data->lmi_map, get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);

      t8_forest_iterate_replace (
        new_forest, forest,
        [] (auto *forest_old, auto *forest_new, auto which_tree, const auto tree_class, const auto *scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, auto first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      // Clean up old forest and user data
      auto *old_user_data = get_user_data ();
      delete old_user_data->lmi_map;
      sc_array_destroy (old_user_data->lmi_idx);
      T8_FREE (old_user_data);
      t8_forest_unref (&forest);

      d_map.erase_all ();
      td_set.erase_all ();
      refinement_set.erase_all ();
      coarsening_set.erase_all ();

      forest = new_forest;
    }
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
