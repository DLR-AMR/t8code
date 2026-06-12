#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/mst.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelindex_set.hxx"
#include "t8_mra/num/dg_basis.hxx"
#include "t8_mra/num/mat.hxx"
#include "t8_cmesh.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"

#include <algorithm>
#include <array>
#include <vector>

namespace t8_mra
{

/**
 * @brief Forward declaration of multiscale_data
 *
 * Stores mask coefficients for MST operations.
 * Should be specialized for each element type.
 */
template <t8_eclass TShape>
struct multiscale_data
{
  static constexpr unsigned short DIM = 0;
  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

/**
 * @brief Forward declaration of multiscale (primary template)
 *
 * This is the primary template that will be specialized for
 * different element types (triangle, quad, etc.)
 */
template <t8_eclass TShape, int U, int P>
class multiscale;

/**
 * @brief Multiscale analysis base class
 *
 * This class template provides a complete MRA implementation that works
 * for both triangular and cartesian elements. It contains all common
 * functionality including:
 *   - Multiscale transformations (forward/inverse)
 *   - Thresholding and adaptation
 *   - Forest management
 *   - Data structure management
 *
 * Element-specific behavior is controlled via policy classes and
 * virtual functions that can be overridden in derived classes.
 *
 * @tparam TShape Element shape (T8_ECLASS_TRIANGLE, T8_ECLASS_QUAD, etc.)
 * @tparam U Number of solution components
 * @tparam P Polynomial order (number of nodes per direction)
 */
template <t8_eclass TShape, unsigned short U, unsigned short P>
class multiscale_base: public multiscale_data<TShape> {
 public:
  using element_t = element_data<TShape, U, P>;
  using levelmultiindex = t8_mra::levelmultiindex<TShape>;
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;
  using MST = mst<element_t>;

  static constexpr auto Shape = TShape;
  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;
  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

  // Bring mask coefficients into scope
  using multiscale_data<TShape>::mask_coefficients;
  using multiscale_data<TShape>::inverse_mask_coefficients;

  //=============================================================================
  // Member Variables (Common to all element types)
  //=============================================================================

  /// Maximum refinement level
  unsigned int maximum_level;

  /// Scaling factors for each solution component (set by criteria via prepare)
  std::array<double, U> c_scaling;

  /// DG basis for projection
  t8_mra::dg_basis<element_t> basis;

  /// Detail coefficient storage
  levelindex_map<levelmultiindex, element_t> d_map;

  /// Set of significant details (thresholding)
  levelindex_set<levelmultiindex> td_set;

  /// Set of elements marked for refinement
  levelindex_set<levelmultiindex> refinement_set;

  /// Set of elements marked for coarsening
  levelindex_set<levelmultiindex> coarsening_set;

  /// t8code forest
  t8_forest_t forest;

  /// MPI communicator
  sc_MPI_Comm comm;

  /// Default quadrature accuracy, derived from the polynomial order: exact
  /// for products of two basis functions (degree 2(P-1)), with margin for
  /// the projection of general data.
  static constexpr int default_dunavant_rule = 2 * P;       // rule number == polynomial exactness
  static constexpr int default_num_quad_points_1d = P + 1;  // n points exact to degree 2n-1

 public:
  //=============================================================================
  // Constructors
  //=============================================================================

  /**
   * @brief Constructor for cartesian elements (QUAD, LINE, HEX)
   *
   * @param _max_level Maximum refinement level
   * @param _comm MPI communicator
   */
  multiscale_base (int _max_level, sc_MPI_Comm _comm)
    requires is_cartesian<TShape>
    : maximum_level (_max_level), basis (default_num_quad_points_1d, P_DIM), d_map (maximum_level),
      td_set (maximum_level), refinement_set (maximum_level), coarsening_set (maximum_level), comm (_comm)
  {
    c_scaling.fill (1.0);
  }

  /**
   * @brief Constructor for triangular elements
   *
   * @param _max_level Maximum refinement level
   * @param _comm MPI communicator
   */
  multiscale_base (int _max_level, sc_MPI_Comm _comm)
    requires (TShape == T8_ECLASS_TRIANGLE)
    : maximum_level (_max_level), basis (t8_mra::dunavant_order_num (default_dunavant_rule), default_dunavant_rule),
      d_map (maximum_level), td_set (maximum_level), refinement_set (maximum_level), coarsening_set (maximum_level),
      comm (_comm)
  {
    c_scaling.fill (1.0);
  }

  virtual ~multiscale_base () = default;

 public:
  //=============================================================================
  // Forest and Data Access
  //=============================================================================

  /**
   * @brief Get the t8code forest
   */
  t8_forest_t
  get_forest ()
  {
    return forest;
  }

  /**
   * @brief Get user data attached to the forest
   */
  t8_mra::forest_data<element_t> *
  get_user_data ()
  {
    return reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest));
  }

  /**
   * @brief Get the level-multiindex map
   */
  t8_mra::levelindex_map<levelmultiindex, element_t> *
  get_lmi_map ()
  {
    return get_user_data ()->lmi_map;
  }

  //=============================================================================
  // Multiscale Transformation
  //=============================================================================

  /**
   * @brief Forward multiscale transformation (restriction: fine -> coarse)
   *
   * Computes parent coefficients and detail coefficients using the
   * MST implementation.
   *
   * @param l_min Minimum refinement level
   * @param l_max Maximum refinement level
   */
  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    MST::forward_transformation (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  /**
   * @brief Inverse multiscale transformation (prolongation: coarse -> fine)
   *
   * Reconstructs children from parent and detail coefficients using the
   * MST implementation.
   *
   * @param l_min Minimum refinement level
   * @param l_max Maximum refinement level
   */
  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    MST::inverse_transformation (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  /**
   * @brief Compute details of leaf families without modifying the grid data
   *
   * Non-destructive counterpart of multiscale_transformation: fills d_map at
   * the parent levels of complete leaf families in (l_min, l_max], while
   * lmi_map keeps the single-scale leaf representation. Basis for the
   * refinement criterion (thresholding / Harten's prediction on details).
   *
   * @param l_min Minimum refinement level
   * @param l_max Maximum refinement level
   */
  void
  compute_leaf_details (unsigned int l_min, unsigned int l_max)
  {
    MST::leaf_details (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  //=============================================================================
  // Thresholding (Element-specific via virtual function)
  //=============================================================================

  /**
   * @brief Compute local detail norm for a given element
   *
   * This function must be implemented by derived classes as the
   * detail norm computation may be element-specific.
   *
   * @param lmi Level multi-index
   * @return Array of detail norms (one per solution component)
   */
  virtual std::array<double, U_DIM>
  local_detail_norm (const levelmultiindex &lmi) = 0;

  /**
   * @brief Maximum detail norm over all components, scaled by c_scaling
   *
   * Common building block for detail-based adaptation criteria.
   *
   * @param lmi Level multi-index
   * @return max_u ||d_u|| / c_scaling_u
   */
  double
  scaled_detail_norm (const levelmultiindex &lmi)
  {
    auto detail_norm = local_detail_norm (lmi);
    for (auto u = 0u; u < U_DIM; ++u)
      detail_norm[u] /= c_scaling[u];

    return *std::max_element (detail_norm.begin (), detail_norm.end ());
  }

  /**
   * @brief Compute local threshold value for an element
   *
   * Implements uniform subdivision thresholding (Veli eq. 2.44)
   *
   * @param lmi Level multi-index
   * @param gamma Expected order of convergence (criterion parameter)
   * @return Local threshold value
   */
  double
  local_threshold_value (const levelmultiindex &lmi, int gamma)
  {
    const auto vol = d_map.get (lmi).vol;

    const auto level_diff = maximum_level - lmi.level ();
    const auto h_lambda = std::sqrt (vol);
    const auto h_max_level = std::pow (vol / std::pow (levelmultiindex::NUM_CHILDREN, level_diff), (gamma + 1.0) / 2.0);

    return h_max_level / h_lambda;
  }

  /**
   * @brief Compute threshold scaling factor (eq. 2.39)
   *
   * Returns scaling factors for each solution component based on
   * mean values over the domain.
   *
   * @return std::array<double, U_DIM> Scaling factor for each component
   */
  std::array<double, U_DIM>
  threshold_scaling_factor ()
  {
    std::array<double, U_DIM> res = {};

    auto current_idx = 0u;
    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    for (auto tree_idx = 0u; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

        const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), current_idx);
        const auto vol = t8_forest_element_volume (forest, tree_idx, element);

        // Compute mean value for each component
        for (auto u = 0u; u < U_DIM; ++u) {
          // Mean value is approximately the first DG coefficient (constant mode)
          // times the scaling function value at the element center
          const auto mean_val = get_lmi_map ()->get (lmi).u_coeffs[element_t::dg_idx (u, 0)];
          res[u] += std::abs (mean_val) * vol;
        }
      }
    }

    // The scaling is a domain integral (eq. 2.39): sum the per-rank
    // partials so all ranks threshold consistently.
    std::array<double, U_DIM> global_res = {};
    sc_MPI_Allreduce (res.data (), global_res.data (), U_DIM, sc_MPI_DOUBLE, sc_MPI_SUM, comm);

    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = std::max (1.0, global_res[u]);

    return res;
  }

  //=============================================================================
  // Projection (Element-specific, must be implemented by derived classes)
  //=============================================================================

  /**
   * @brief Project a function onto the DG basis for an element
   *
   * This function must be implemented by derived classes as projection
   * is element-specific (different quadrature, coordinate mappings, etc.)
   *
   * NOTE: This is a template method in derived classes and cannot be virtual.
   * Templates cannot be virtual in C++.
   *
   * Derived classes should implement:
   *   template <typename Func>
   *   void project_impl(std::vector<double> &dg_coeffs, int tree_idx,
   *                     const t8_element_t *element, Func &&func)
   */

  //=============================================================================
  // Post-Adaptation Hook
  //=============================================================================

  /**
   * @brief Post-adaptation hook (for element-specific operations)
   *
   * Override this in derived classes to perform element-specific operations
   * after forest adaptation (e.g., update vertex orders for triangles).
   *
   * Default implementation does nothing (for cartesian elements).
   */
  virtual void
  post_adaptation_hook ()
  {
    // Default: no-op
  }

  //=============================================================================
  // Cleanup
  //=============================================================================

  /**
   * @brief Clean up all data structures
   */
  void
  cleanup ()
  {
    d_map.erase_all ();
    td_set.erase_all ();
    refinement_set.erase_all ();
    coarsening_set.erase_all ();
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
