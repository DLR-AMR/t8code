#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/mst.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelindex_set.hxx"
#include "t8_mra/num/dg_basis.hxx"
#include "t8_mra/num/mat.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_forest/t8_forest_iterate.h"

#include <algorithm>
#include <array>
#include <optional>
#include <vector>

namespace t8_mra
{

/// Per-shape mask coefficients for the MST.
template <t8_eclass TShape>
struct multiscale_data
{
  static constexpr unsigned short DIM = 0;
  std::vector<t8_mra::mat> mask_coefficients;
};

/// Primary template; specialized per shape in shapes/.
template <t8_eclass TShape, int U, int P>
class multiscale;

/**
 * @brief Common MRA functionality for all element shapes: the multiscale
 * transformations, thresholding, projection, and forest/data management.
 * Shape-specific behaviour is supplied by the derived multiscale<> through
 * CRTP hooks (local_detail_norm, evaluate, evaluate_gradient, project_impl,
 * post_adaptation_hook), resolved statically via derived().
 *
 * @tparam Derived the multiscale<> specialization (CRTP)
 * @tparam TShape element shape, @tparam U solution components, @tparam P order
 */
template <typename Derived, t8_eclass TShape, unsigned short U, unsigned short P>
class multiscale_base: public multiscale_data<TShape> {
 public:
  Derived &
  derived ()
  {
    return static_cast<Derived &> (*this);
  }

  const Derived &
  derived () const
  {
    return static_cast<const Derived &> (*this);
  }

  using element_t = element_data<TShape, U, P>;
  using detail_t = detail_data<TShape, U, P>;
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
  levelindex_map<levelmultiindex, detail_t> d_map;

  /// Set of significant details (thresholding)
  levelindex_set<levelmultiindex> td_set;

  /// Set of elements marked for refinement
  levelindex_set<levelmultiindex> refinement_set;

  /// Set of elements marked for coarsening
  levelindex_set<levelmultiindex> coarsening_set;

  /// Ghost element data, keyed by lmi. Read-only snapshot filled by
  /// ghost_exchange; invalidated by every adaptation/repartition. Kept
  /// separate from lmi_map, where ghost entries would masquerade as local
  /// leaves in the family analysis.
  levelindex_map<levelmultiindex, element_t> ghost_map;

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
    : maximum_level (_max_level), basis (default_num_quad_points_1d), d_map (maximum_level), td_set (maximum_level),
      refinement_set (maximum_level), coarsening_set (maximum_level), ghost_map (maximum_level), comm (_comm)
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
    : maximum_level (_max_level), basis (default_dunavant_rule), d_map (maximum_level), td_set (maximum_level),
      refinement_set (maximum_level), coarsening_set (maximum_level), ghost_map (maximum_level), comm (_comm)
  {
    c_scaling.fill (1.0);
  }

  ~multiscale_base () = default;

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
  // Local leaf iteration / data construction
  //=============================================================================

  /// Iterate the local leaves in SFC order. f receives (local tree index,
  /// element, forest-local leaf index, global tree id); the leaf index matches
  /// the lmi_idx array layout.
  ///
  /// f: (t8_locidx_t tree_idx, const t8_element_t *element,
  ///     unsigned int local_idx, t8_gloidx_t global_tree) -> void
  template <typename F>
  void
  for_each_local_leaf (F &&f)
  {
    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto local_idx = 0u;
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elems = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);
      const auto global_tree = t8_forest_global_tree_id (forest, tree_idx);
      for (t8_locidx_t ele_idx = 0; ele_idx < num_elems; ++ele_idx, ++local_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        f (tree_idx, element, local_idx, global_tree);
      }
    }
  }

  /// Build lmi_map + lmi_idx for the committed forest and attach as user data.
  /// projector supplies the per-leaf data.
  ///
  /// projector: (int tree_idx, const t8_element_t *element) -> element_t
  template <typename Projector>
  void
  build_lmi_map (const t8_scheme *scheme, Projector &&projector)
  {
    T8_ASSERT (t8_forest_is_committed (forest));

    auto *user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    const auto num_local = t8_forest_get_local_num_leaf_elements (forest);
    const auto num_ghost = t8_forest_get_num_ghosts (forest);

    user_data->lmi_map = new levelindex_map<levelmultiindex, element_t> (maximum_level);
    user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local + num_ghost);
    user_data->mra_instance = &derived ();
    t8_forest_set_user_data (forest, user_data);

    for_each_local_leaf (
      [&] (t8_locidx_t tree_idx, const t8_element_t *element, unsigned int local_idx, t8_gloidx_t global_tree) {
        const auto lmi = levelmultiindex (global_tree, element, scheme);
        user_data->lmi_map->insert (lmi, projector (tree_idx, element));
        t8_mra::set_lmi_forest_data (user_data, local_idx, lmi);
      });
  }

  //=============================================================================
  // Multiscale Transformation
  //=============================================================================

  /**
   * @brief Forward multiscale transformation (analysis), non-destructive.
   *
   * Fills d_map with parent coeffs + details of complete families in
   * (l_min, l_max]; lmi_map keeps its single-scale leaves. This is the
   * analysis coarsen and refine run to obtain the details their criteria
   * act on.
   */
  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    MST::multiscale_transformation (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  /**
   * @brief Inverse multiscale transformation (reconstruction: coarse -> fine).
   *
   * Reconstructs children from parent + details. Refinement's realization step
   * (children from zero details = exact subdivision).
   */
  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    MST::inverse_multiscale_transformation (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  /**
   * @brief Full forward multiscale transformation (restriction), destructive.
   *
   * Collapses the single-scale representation to l_min (coarsest scaling
   * coeffs + details). Inverse of inverse_multiscale_transformation; used for
   * the MST round-trip property.
   */
  void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max)
  {
    MST::multiscale_decomposition (l_min, l_max, get_lmi_map (), d_map, mask_coefficients);
  }

  //=============================================================================
  // Thresholding
  //=============================================================================

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
    auto detail_norm = derived ().local_detail_norm (lmi);
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

    for_each_local_leaf ([&] (t8_locidx_t, const t8_element_t *, unsigned int local_idx, t8_gloidx_t) {
      const auto lmi = t8_mra::get_lmi_from_forest_data (get_user_data (), local_idx);
      const auto &data = get_lmi_map ()->get (lmi);
      const auto mean = mean_val (data);

      for (auto u = 0u; u < U_DIM; ++u)
        res[u] += std::abs (mean[u]) * data.vol;
    });

    // The scaling is a domain integral (eq. 2.39): sum the per-rank
    // partials so all ranks threshold consistently.
    std::array<double, U_DIM> global_res = {};
    sc_MPI_Allreduce (res.data (), global_res.data (), U_DIM, sc_MPI_DOUBLE, sc_MPI_SUM, comm);

    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = std::max (1.0, global_res[u]);

    return res;
  }

  //=============================================================================
  // Evaluation
  //=============================================================================

  /**
   * @brief Evaluate the solution at a reference-cell point
   *
   * x_ref is [0,1]^DIM (cartesian) or (tau1, tau2) on the reference triangle;
   * the triangle -> basis {lambda0, lambda1} conversion is handled here.
   */
  std::array<double, U_DIM>
  evaluate_reference (const element_t &data, const std::array<double, DIM> &x_ref)
  {
    std::array<double, DIM> x_basis = x_ref;
    if constexpr (Shape == T8_ECLASS_TRIANGLE)
      x_basis = { 1.0 - x_ref[0] - x_ref[1], x_ref[0] };  // (tau1, tau2) -> (lambda0, lambda1)

    const std::vector<double> x (x_basis.begin (), x_basis.end ());
    const auto phi = basis.basis_value (x);
    const auto scaling = t8_mra::basis<Shape, P_DIM>::normalization (data.vol);

    std::array<double, U_DIM> res = {};
    for (auto u = 0u; u < U_DIM; ++u)
      for (auto i = 0u; i < DOF; ++i)
        res[u] += data.u_coeffs[element_t::dg_idx (u, i)] * scaling * phi[i];

    return res;
  }

  /// A point-location query for t8_forest_search, filled with the value at the
  /// owning leaf. Trivially copyable to live in the search's sc_array.
  struct point_query
  {
    double point[3];
    double tolerance;
    int found;
    std::array<double, U_DIM> value;
  };

  /// Search criterion: descend everywhere, the queries do the pruning.
  static int
  search_descend_fn (t8_forest_t, const t8_locidx_t, const t8_element_t *, const int, const t8_element_array_t *,
                     const t8_locidx_t)
  {
    return 1;
  }

  /// Per-element query: stay active while the point is inside (so the search
  /// recurses to the owner), evaluate at the owning leaf. MRA via forest user
  /// data.
  static void
  search_point_fn (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                   const t8_element_array_t *, const t8_locidx_t, sc_array_t *queries, sc_array_t *query_indices,
                   int *query_matches, const size_t num_active_queries)
  {
    auto *user_data = reinterpret_cast<forest_data<element_t> *> (t8_forest_get_user_data (forest));
    auto *mra = static_cast<Derived *> (user_data->mra_instance);
    const auto *scheme = t8_forest_get_scheme (forest);

    for (size_t i = 0; i < num_active_queries; ++i) {
      const size_t query_idx = *static_cast<size_t *> (sc_array_index (query_indices, i));
      auto *query = static_cast<point_query *> (sc_array_index (queries, query_idx));

      int inside = 0;
      t8_forest_element_points_inside (forest, ltreeid, element, query->point, 1, &inside, query->tolerance);
      query_matches[i] = inside;

      if (inside && is_leaf && !query->found) {
        const auto gtree = t8_forest_global_tree_id (forest, ltreeid);
        const auto lmi = levelmultiindex (gtree, element, scheme);

        std::array<double, DIM> x;
        for (auto d = 0u; d < DIM; ++d)
          x[d] = query->point[d];

        query->value = mra->evaluate (ltreeid, element, mra->get_lmi_map ()->get (lmi), x);
        query->found = 1;
      }
    }
  }

  /**
   * @brief Evaluate the solution at an arbitrary physical point of the domain
   *
   * Locates the owning local leaf via t8_forest_search; nullopt if no local leaf
   * owns the point (outside the domain, or on another rank).
   */
  std::optional<std::array<double, U_DIM>>
  evaluate_point (const std::array<double, DIM> &x, double tolerance = 1e-8)
  {
    point_query query = {};
    for (auto d = 0u; d < DIM; ++d)
      query.point[d] = x[d];
    query.tolerance = tolerance;

    sc_array_t *queries = sc_array_new_count (sizeof (point_query), 1);
    *static_cast<point_query *> (sc_array_index (queries, 0)) = query;

    t8_forest_search (forest, search_descend_fn, search_point_fn, queries);

    const point_query result = *static_cast<point_query *> (sc_array_index (queries, 0));
    sc_array_destroy (queries);

    if (result.found)
      return result.value;
    return std::nullopt;
  }

  /// Cell average per component.
  std::array<double, U_DIM>
  mean_val (const element_t &data)
  {
    const double scale = ((Shape == T8_ECLASS_TRIANGLE) ? std::sqrt (data.vol) : data.vol) / data.vol;

    std::array<double, U_DIM> mean;
    for (auto u = 0u; u < U_DIM; ++u)
      mean[u] = data.u_coeffs[element_t::dg_idx (u, 0)] * scale;

    return mean;
  }

  /// Cell average of the leaf lmi.
  std::array<double, U_DIM>
  mean_val (const levelmultiindex &lmi)
  {
    return mean_val (get_lmi_map ()->get (lmi));
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
   * A template method in the derived class
   *
   * Derived classes should implement:
   *   template <typename Func>
   *   void project_impl(std::span<double> dg_coeffs, int tree_idx,
   *                     const t8_element_t *element, Func &&func)
   */

  //=============================================================================
  // Post-Adaptation Hook
  //=============================================================================

  /**
   * @brief Post-adaptation hook
   *
   * Runs after every forest adaptation. The derived shape shadows this to do
   * element-specific work (triangle updates vertex orders); the default is a
   * no-op (cartesian).
   */
  void
  post_adaptation_hook ()
  {
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
    ghost_map.erase_all ();

    if (forest != nullptr) {
      if (auto *user_data = get_user_data ()) {
        delete user_data->lmi_map;
        if (user_data->lmi_idx)
          sc_array_destroy (user_data->lmi_idx);
        T8_FREE (user_data);
      }
      t8_forest_unref (&forest);
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
