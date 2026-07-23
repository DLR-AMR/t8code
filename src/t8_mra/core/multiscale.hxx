#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/forest_backend.hxx"
#include "t8_mra/core/mst.hxx"
#include "t8_mra/core/adapt/coarsen.hxx"
#include "t8_mra/core/adapt/refine.hxx"
#include "t8_mra/core/adapt/balance.hxx"
#include "t8_mra/dg/cartesian.hxx"
#include "t8_mra/dg/triangle.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelindex_set.hxx"
#include "t8_mra/criteria/coarsening_criterion.hxx"
#include "t8_mra/criteria/refinement_criterion.hxx"
#include "t8_mra/num/cell_geometry.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/nodal_to_modal.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_iterate.h"

#include <algorithm>
#include <array>
#include <optional>
#include <span>

namespace t8_mra
{

/**
 * @brief Multiresolution analysis on an adaptive t8code forest.
 *
 * Composition of three collaborators:
 *   - grid:      the t8code + MPI side (forest, lmi_map, ghost, adapt/partition)
 *   - transform: the multiscale (two-scale) transform on the maps
 *   - dg:        the per-shape DG numerics (projection, evaluation, geometry)
 *
 * The adaptation operations (coarsen/refine/balance) live as free functions in
 * core/adapt/ and are forwarded here.
 *
 * @tparam TShape element shape, @tparam U components, @tparam P order
 */
template <t8_eclass TShape, int U, int P>
class multiscale {
 public:
  using element_t = element_data<TShape, U, P>;
  using detail_t = detail_data<TShape, U, P>;
  using levelmultiindex = t8_mra::levelmultiindex<TShape>;
  using index_set = ankerl::unordered_dense::set<levelmultiindex>;
  using geometry_t = cell_geometry<TShape, P>;
  using dg_t = dg<TShape, U, P>;
  using MST = mst<element_t>;

  static constexpr auto Shape = TShape;
  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;
  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

  //=============================================================================
  // Collaborators and state
  //=============================================================================

  /// t8code + MPI side
  forest_backend<TShape, U, P> grid;

  /// Two-scale transform (owns the mask coefficients)
  MST transform;

  /// Per-shape DG numerics
  dg_t discretization;

  /// Scaling factors per component (set by criteria via prepare)
  std::array<double, U> c_scaling;

  /// Detail coefficient storage
  levelindex_map<levelmultiindex, detail_t> d_map;

  /// Significant details
  levelindex_set<levelmultiindex> td_set;

  /// Elements marked for refinement
  levelindex_set<levelmultiindex> refinement_set;

  /// Elements marked for coarsening
  levelindex_set<levelmultiindex> coarsening_set;

  multiscale (int max_level, sc_MPI_Comm comm)
    : grid (max_level, comm), d_map (grid.maximum_level), td_set (grid.maximum_level),
      refinement_set (grid.maximum_level), coarsening_set (grid.maximum_level)
  {
    c_scaling.fill (1.0);
    grid.bind (this, [this] () { post_adapt (); });
  }

  //=============================================================================
  // Accessors (forward to grid)
  //=============================================================================

  t8_forest_t
  get_forest ()
  {
    return grid.get_forest ();
  }

  sc_MPI_Comm
  get_comm () const
  {
    return grid.comm;
  }

  t8_mra::forest_data<element_t> *
  get_user_data ()
  {
    return grid.get_user_data ();
  }

  levelindex_map<levelmultiindex, element_t> *
  get_lmi_map ()
  {
    return grid.get_lmi_map ();
  }

  unsigned int
  maximum_level () const
  {
    return grid.maximum_level;
  }

  template <typename F>
  void
  for_each_local_leaf (F &&f)
  {
    grid.for_each_local_leaf (std::forward<F> (f));
  }

  void
  ghost_exchange ()
  {
    grid.ghost_exchange ();
  }

  void
  repartition ()
  {
    grid.repartition ();
  }

  //=============================================================================
  // Multiscale transform (forward to transform)
  //=============================================================================

  void
  multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    transform.multiscale_transformation (l_min, l_max, get_lmi_map (), d_map);
  }

  void
  inverse_multiscale_transformation (unsigned int l_min, unsigned int l_max)
  {
    transform.inverse_multiscale_transformation (l_min, l_max, get_lmi_map (), d_map);
  }

  void
  multiscale_decomposition (unsigned int l_min, unsigned int l_max)
  {
    transform.multiscale_decomposition (l_min, l_max, get_lmi_map (), d_map);
  }

  //=============================================================================
  // Thresholding
  //=============================================================================

  /** @brief max_u ||d_u|| / c_scaling_u. */
  double
  scaled_detail_norm (const levelmultiindex &lmi)
  {
    auto detail_norm = transform.detail_norm (d_map.get (lmi));
    for (auto u = 0u; u < U_DIM; ++u)
      detail_norm[u] /= c_scaling[u];

    return *std::max_element (detail_norm.begin (), detail_norm.end ());
  }

  /** @brief Level-dependent threshold (Veli eq. 2.44). */
  double
  local_threshold_value (const levelmultiindex &lmi, int gamma)
  {
    const auto vol = d_map.get (lmi).vol;

    const auto level_diff = grid.maximum_level - lmi.level ();
    const auto h_lambda = std::sqrt (vol);
    const auto h_max_level = std::pow (vol / std::pow (levelmultiindex::NUM_CHILDREN, level_diff), (gamma + 1.0) / 2.0);

    return h_max_level / h_lambda;
  }

  /** @brief Per-component domain-integral scaling (eq. 2.39), reduced over ranks. */
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

    std::array<double, U_DIM> global_res = {};
    sc_MPI_Allreduce (res.data (), global_res.data (), U_DIM, sc_MPI_DOUBLE, sc_MPI_SUM, grid.comm);

    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = std::max (1.0, global_res[u]);

    return res;
  }

  //=============================================================================
  // Evaluation
  //=============================================================================

  /** @brief Solution value per component at a reference-cell point. */
  std::array<double, U_DIM>
  evaluate_reference (const element_t &data, const std::array<double, DIM> &x_ref)
  {
    std::array<double, U_DIM> res = {};
    for (auto u = 0u; u < U_DIM; ++u)
      res[u] = geometry_t::reference_value (
        std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), x_ref, data.vol);

    return res;
  }

  /** @brief Cell average per component. */
  std::array<double, U_DIM>
  mean_val (const element_t &data)
  {
    std::array<double, U_DIM> mean;
    for (auto u = 0u; u < U_DIM; ++u)
      mean[u] = t8_mra::cell_mean<TShape, P> (
        std::span<const double> (&data.u_coeffs[element_t::dg_idx (u, 0)], DOF), data.vol);

    return mean;
  }

  std::array<double, U_DIM>
  mean_val (const levelmultiindex &lmi)
  {
    return mean_val (get_lmi_map ()->get (lmi));
  }

  /** @brief Solution value per component at a physical point of the given leaf. */
  std::array<double, U_DIM>
  evaluate (int tree_idx, const t8_element_t *element, const element_t &data, const std::array<double, DIM> &x_phys)
  {
    double corners[T8_ECLASS_MAX_CORNERS][3];
    grid.element_corner_coords (tree_idx, element, corners);
    const auto geom = discretization.geometry (corners, data.vol, data.order);
    return discretization.evaluate (geom, data, x_phys);
  }

  /** @brief Solution gradient per component at a physical point of the given leaf. */
  std::array<std::array<double, DIM>, U_DIM>
  evaluate_gradient (int tree_idx, const t8_element_t *element, const element_t &data,
                     const std::array<double, DIM> &x_phys)
  {
    double corners[T8_ECLASS_MAX_CORNERS][3];
    grid.element_corner_coords (tree_idx, element, corners);
    const auto geom = discretization.geometry (corners, data.vol, data.order);
    return discretization.evaluate_gradient (geom, data, x_phys);
  }

  /// Point-location query for t8_forest_search, filled with the owning leaf value.
  struct point_query
  {
    double point[3];
    double tolerance;
    int found;
    std::array<double, U_DIM> value;
  };

  static int
  search_descend_fn (t8_forest_t, const t8_locidx_t, const t8_element_t *, const int, const t8_element_array_t *,
                     const t8_locidx_t)
  {
    return 1;
  }

  static void
  search_point_fn (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                   const t8_element_array_t *, const t8_locidx_t, sc_array_t *queries, sc_array_t *query_indices,
                   int *query_matches, const size_t num_active_queries)
  {
    auto *user_data = reinterpret_cast<forest_data<element_t> *> (t8_forest_get_user_data (forest));
    auto *mra = static_cast<multiscale *> (user_data->mra_instance);
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

  /** @brief Solution value at a physical point; nullopt if no local leaf owns it. */
  std::optional<std::array<double, U_DIM>>
  evaluate_point (const std::array<double, DIM> &x, double tolerance = 1e-8)
  {
    point_query query = {};
    for (auto d = 0u; d < DIM; ++d)
      query.point[d] = x[d];
    query.tolerance = tolerance;

    sc_array_t *queries = sc_array_new_count (sizeof (point_query), 1);
    *static_cast<point_query *> (sc_array_index (queries, 0)) = query;

    t8_forest_search (grid.get_forest (), search_descend_fn, search_point_fn, queries);

    const point_query result = *static_cast<point_query *> (sc_array_index (queries, 0));
    sc_array_destroy (queries);

    if (result.found)
      return result.value;
    return std::nullopt;
  }

  //=============================================================================
  // Projection / initialization
  //=============================================================================

  /** @brief Project func onto a single leaf. */
  template <typename Func>
  element_t
  project_leaf (int tree_idx, const t8_element_t *element, Func &&func)
  {
    element_t data;
    data.vol = grid.element_volume (tree_idx, element);

    const auto *scheme = t8_forest_get_scheme (grid.get_forest ());
    data.order = levelmultiindex::point_order_at_level (element, scheme);

    double corners[T8_ECLASS_MAX_CORNERS][3];
    grid.element_corner_coords (tree_idx, element, corners);
    const auto geom = discretization.geometry (corners, data.vol, data.order);
    discretization.project (data.u_coeffs, geom, func);

    return data;
  }

  /** @brief Project func onto a uniform forest of the given level. */
  template <typename Func>
  void
  initialize_data (t8_cmesh_t mesh, const t8_scheme *scheme, int level, Func &&func)
  {
    grid.forest = t8_forest_new_uniform (mesh, scheme, level, 0, grid.comm);
    grid.build ([&] (int tree_idx, const t8_element_t *element) { return project_leaf (tree_idx, element, func); });
  }

  /** @brief Load per-cell nodal DG values onto an existing forest as modal coeffs. */
  template <typename CellNodalValues>
  void
  initialize_data_nodal (t8_forest_t forest, const std::array<std::array<double, DIM>, DOF> &nodes,
                         CellNodalValues &&cell_nodal_values)
  {
    const nodal_to_modal<TShape, U, P> to_modal (nodes);

    t8_forest_ref (forest);
    grid.forest = forest;

    grid.build ([&] (int tree_idx, const t8_element_t *element) {
      element_t data;
      data.vol = grid.element_volume (tree_idx, element);
      const auto *scheme = t8_forest_get_scheme (grid.get_forest ());
      data.order = levelmultiindex::point_order_at_level (element, scheme);
      const auto nodal = cell_nodal_values (tree_idx, element);
      to_modal (std::span<const double> (nodal.data (), nodal.size ()), data.u_coeffs);

      return data;
    });
  }

  /** @brief Reconstruct per-cell nodal DG values from the current forest. */
  template <typename WriteCellNodalValues>
  void
  export_data_nodal (const std::array<std::array<double, DIM>, DOF> &nodes, WriteCellNodalValues &&write_cell_nodal_values)
  {
    const modal_to_nodal<TShape, U, P> to_nodal (nodes);
    const auto *scheme = t8_forest_get_scheme (grid.get_forest ());
    auto *lmi_map = get_lmi_map ();

    for_each_local_leaf ([&] (t8_locidx_t tree_idx, const t8_element_t *element, unsigned int, t8_gloidx_t global_tree) {
      const levelmultiindex lmi (global_tree, element, scheme);
      const auto *data = lmi_map->find (lmi);
      const auto nodal = to_nodal (std::span<const double> (data->u_coeffs.data (), data->u_coeffs.size ()));
      write_cell_nodal_values (tree_idx, element, std::span<const double> (nodal.data (), nodal.size ()));
    });
  }

  //=============================================================================
  // Adaptation (forward to adapt::)
  //=============================================================================

  template <typename Criterion = hard_thresholding>
    requires coarsening_criterion<Criterion, multiscale>
  void
  coarsen (int min_level, int max_level, Criterion criterion = {})
  {
    adapt::coarsen (*this, min_level, max_level, criterion);
  }

  template <typename Criterion = harten_prediction>
    requires refinement_criterion<Criterion, multiscale>
  void
  refine (int min_level, int max_level, Criterion criterion = {})
  {
    adapt::refine (*this, min_level, max_level, criterion);
  }

  void
  balance ()
  {
    adapt::balance (*this);
  }

  template <typename Func, typename Criterion = hard_thresholding>
    requires coarsening_criterion<Criterion, multiscale>
  void
  initialize_data_adaptive (t8_cmesh_t mesh, const t8_scheme *scheme, int max_level, Func &&func, Criterion criterion = {})
  {
    adapt::initialize_data_adaptive (*this, mesh, scheme, max_level, func, criterion);
  }

  //=============================================================================
  // t8code adaptation callbacks
  //=============================================================================

  int
  coarsening_callback (t8_forest_t, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t, t8_locidx_t,
                       const t8_scheme_c *scheme, int is_family, int, t8_element_t *elements[])
  {
    if (!is_family)
      return 0;

    const auto gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
    const auto lmi = levelmultiindex (gtreeid, elements[0], scheme);

    return coarsening_set.contains (lmi) ? -1 : 0;
  }

  int
  refinement_callback (t8_forest_t, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t, t8_locidx_t,
                       const t8_scheme_c *scheme, int, int, t8_element_t *elements[])
  {
    const auto gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
    const auto lmi = levelmultiindex (gtreeid, elements[0], scheme);

    return refinement_set.contains (lmi) ? 1 : 0;
  }

  static int
  static_coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    auto *user_data = reinterpret_cast<forest_data<element_t> *> (t8_forest_get_user_data (forest_from));
    return static_cast<multiscale *> (user_data->mra_instance)
      ->coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                             num_elements, elements);
  }

  static int
  static_refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    auto *user_data = reinterpret_cast<forest_data<element_t> *> (t8_forest_get_user_data (forest_from));
    return static_cast<multiscale *> (user_data->mra_instance)
      ->refinement_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family,
                             num_elements, elements);
  }

  //=============================================================================
  // Post-adaptation hook + cleanup
  //=============================================================================

  /** @brief Refresh the per-leaf vertex order (triangle Bey type); no-op values for cartesian. */
  void
  post_adapt ()
  {
    if (grid.get_forest () == nullptr)
      return;

    auto *user_data = grid.get_user_data ();
    const auto *scheme = t8_forest_get_scheme (grid.get_forest ());

    grid.for_each_local_leaf ([&] (t8_locidx_t, const t8_element_t *elem, unsigned int local_idx, t8_gloidx_t) {
      const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_idx);
      if (auto *data = user_data->lmi_map->find (lmi))
        data->order = levelmultiindex::point_order_at_level (elem, scheme);
    });
  }

  void
  cleanup ()
  {
    d_map.erase_all ();
    td_set.erase_all ();
    refinement_set.erase_all ();
    coarsening_set.erase_all ();
    grid.cleanup ();
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
