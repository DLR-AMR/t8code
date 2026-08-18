#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <unordered_map>

#include "sc_mpi.h"

#include "t8.h"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_eclass/t8_eclass.h"
#include "t8_element/t8_element.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_schemes/t8_scheme.hxx"

#include "t8_mra/core/adapt/grading.hxx"
#include "t8_mra/criteria/coarsening_criterion.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra::adapt
{

/**
 * @brief Destructive fine->coarse sweep collapsing non-significant families
 *
 * Thresholded multiscale decomposition: a collapsed parent is a leaf for the
 * next level, so one traversal captures the whole cascade.
 *
 * @return Number of families collapsed on this rank
 */
template <typename TMultiscale, typename TCriterion>
[[nodiscard]] unsigned int
coarsen_sweep (TMultiscale &mra, int min_level, int max_level, TCriterion &criterion)
{
  using levelmultiindex = typename TMultiscale::levelmultiindex;

  auto num_collapsed_children = 0u;
  mra.transform.multiscale_decomposition (
    min_level, max_level, *mra.get_lmi_map (), mra.d_map,
    [&] (const auto &parent) { return criterion.significant (mra, parent); },
    [&] (const auto &child) {
      mra.coarsening_set.insert (child);
      ++num_collapsed_children;
    });

  return num_collapsed_children / levelmultiindex::NUM_CHILDREN;
}

/**
 * @brief Restore the leaves the sweep collapsed across a 2:1 violation
 *
 * Same rounds as balance(), but the children are reconstructed from the details
 * the sweep left in d_map instead of a zero detail, so a restored leaf carries
 * its exact data. Runs on the forest of its own pass, before any repartition can
 * move a cell away from the rank holding its details.
 *
 * @return Number of cells restored on this rank
 */
template <typename TMultiscale>
unsigned int
restore_graded_leaves (TMultiscale &mra)
{
  using levelmultiindex = typename TMultiscale::levelmultiindex;
  using detail_t = typename TMultiscale::detail_t;

  const auto max_level = mra.grid.maximum_level;
  auto num_restored = 0u;

  for (;;) {
    mra.refinement_set.erase_all ();

    const auto num_marked
      = mra.grid.global_num_marks (grade_neighbours (mra, 0, 1u, [] (const auto &lmi) { return lmi.level () >= 2; }));

    if (num_marked == 0)
      break;

    levelindex_map<levelmultiindex, detail_t> marked_details (max_level);
    for (auto l = 0u; l < max_level; ++l)
      for (const auto &lmi : mra.refinement_set[l]) {
        const auto *details = mra.d_map.find (lmi);
        marked_details.insert (lmi, details != nullptr ? *details : detail_t {});
      }

    mra.transform.inverse_multiscale_transformation (0, max_level, *mra.get_lmi_map (), marked_details);
    mra.grid.adapt (TMultiscale::static_refinement_callback, 1);

    num_restored += num_marked;
  }

  mra.refinement_set.erase_all ();

  return num_restored;
}

/**
 * @brief Adaptive coarsening from max_level down to min_level
 *
 * One destructive fine->coarse sweep on the maps per pass; across ranks an outer
 * fixpoint (adapt + repartition make seam families whole) until no rank marks.
 * With graded set, every pass restores the leaves it collapsed across a 2:1
 * violation from its own details, so they carry their exact data instead of a
 * projection. A pass that collapses no more than it restores ends the fixpoint.
 */
template <typename TMultiscale, typename TCriterion = hard_thresholding>
  requires coarsening_criterion<TCriterion, TMultiscale>
void
coarsen (TMultiscale &mra, int min_level, int max_level, TCriterion criterion = {}, bool graded = false)
{
  if constexpr (criterion_has_prepare<TCriterion, TMultiscale>)
    criterion.prepare (mra);

  for (auto pass = 0;; ++pass) {
    clear_state (mra);

    const auto num_leaves = mra.grid.global_num_leaves ();
    const auto num_marked = coarsen_sweep (mra, min_level, max_level, criterion);

    t8_debugf ("MRA coarsen pass %d: %u families marked, %zu leaves remain\n", pass, num_marked,
               mra.get_lmi_map ()->size ());

    if (mra.grid.global_num_marks (num_marked) == 0)
      break;

    mra.grid.adapt (TMultiscale::static_coarsening_callback, 1);

    if (graded)
      restore_graded_leaves (mra);

    mra.grid.repartition ();

    if (graded && mra.grid.global_num_leaves () >= num_leaves)
      break;
  }

  clear_state (mra);
}

/// Per-component max mean magnitude over a level's leaves, reduced across ranks
/// (floored at 1) so every rank normalizes jump detection identically.
template <typename TMultiscale>
[[nodiscard]] auto
global_v_max (TMultiscale &mra, int level)
{
  std::array<double, TMultiscale::U_DIM> local;
  local.fill (1.0);

  for (const auto &[lmi, _] : (*mra.get_lmi_map ())[level]) {
    const auto m = mra.mean_val (lmi);

    for (auto u = 0u; u < TMultiscale::U_DIM; ++u)
      local[u] = std::max (local[u], std::abs (m[u]));
  }

  std::array<double, TMultiscale::U_DIM> global;
  sc_MPI_Allreduce (local.data (), global.data (), TMultiscale::U_DIM, sc_MPI_DOUBLE, sc_MPI_MAX, mra.grid.comm);

  return global;
}

/**
 * @brief Mean-value jump detection on the leaves of one level
 *
 * Marks a family when a face-neighbour mean difference exceeds c_thresh*sqrt(h):
 * smooth data decays as O(h) and falls below, a discontinuity stays O(1). Remote
 * neighbours come from the ghost layer; the result is globalized because coarsen
 * repartitions between passes.
 *
 * @return Parent lmis of the jumping families
 */
template <typename TMultiscale>
[[nodiscard]] auto
detect_jumps (TMultiscale &mra, int level, double c_thresh)
{
  using levelmultiindex = typename TMultiscale::levelmultiindex;

  mra.grid.ghost_exchange ();

  auto *lmi_map = mra.get_lmi_map ();
  auto &ghost_map = mra.grid.ghost_map;
  const auto v_max = global_v_max (mra, level);

  std::unordered_map<size_t, double> face_jump;
  mra.grid.for_each_face_neigh ([&] (const auto &lmi) { return lmi.level () == static_cast<unsigned int> (level); },
                                [&] (const auto &lmi, t8_eclass_t, t8_gloidx_t, t8_element_t *, const auto &neigh_lmi) {
                                  const auto *neigh_data = lmi_map->contains (neigh_lmi)    ? &lmi_map->get (neigh_lmi)
                                                           : ghost_map.contains (neigh_lmi) ? &ghost_map.get (neigh_lmi)
                                                                                            : nullptr;
                                  if (neigh_data == nullptr)
                                    return;

                                  const auto mean_inner = mra.mean_val (lmi_map->get (lmi));
                                  const auto mean_neigh = mra.mean_val (*neigh_data);
                                  auto &diff = face_jump[lmi.index];

                                  for (auto u = 0u; u < TMultiscale::U_DIM; ++u)
                                    diff = std::max (diff, std::abs (mean_inner[u] - mean_neigh[u]) / v_max[u]);
                                });

  typename TMultiscale::index_set jumps;
  for (const auto &[index, diff] : face_jump) {
    const auto lmi = levelmultiindex (index);
    const auto h = std::pow (lmi_map->get (lmi).vol, 1.0 / TMultiscale::DIM);

    if (diff > c_thresh * std::sqrt (h))
      jumps.insert (t8_mra::parent_lmi (lmi));
  }

  mra.grid.globalize (jumps);

  return jumps;
}

/// Coarsening criterion wrapper: families with a detected jump are always
/// significant.
template <typename TCriterion, typename TMultiscale>
struct jump_guarded
{
  TCriterion &criterion;
  const typename TMultiscale::index_set &jumps;

  void
  prepare (TMultiscale &mra)
  {
    if constexpr (criterion_has_prepare<TCriterion, TMultiscale>)
      criterion.prepare (mra);
  }

  bool
  significant (TMultiscale &mra, const typename TMultiscale::levelmultiindex &lmi)
  {
    return jumps.contains (lmi) || criterion.significant (mra, lmi);
  }
};

/**
 * @brief Refine every leaf at the given level and project the initial data
 *
 * Unlike refine(), the children data is projected directly from the initial
 * data (exact up to quadrature), not predicted. 
 *
 * @return Number of leaves refined
 */
template <typename TMultiscale, typename TFunc>
unsigned int
project_onto_children (TMultiscale &mra, int level, TFunc &&func)
{
  clear_state (mra);

  for (const auto &[lmi, _] : (*mra.get_lmi_map ())[level])
    mra.refinement_set.insert (lmi);

  const auto num_marked = mra.refinement_set[level].size ();
  if (mra.grid.global_num_marks (num_marked) == 0)
    return 0;

  mra.grid.adapt (TMultiscale::static_refinement_callback);

  auto *lmi_map = mra.get_lmi_map ();
  auto *user_data = mra.get_user_data ();

  const auto num_local_trees = t8_forest_get_num_local_trees (mra.grid.get_forest ());
  auto current_idx = 0;

  for (auto tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (mra.grid.get_forest (), tree_idx);

    for (auto ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
      const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);

      if (lmi_map->contains (lmi))
        continue;

      const auto *element = t8_forest_get_leaf_element_in_tree (mra.grid.get_forest (), tree_idx, ele_idx);
      lmi_map->insert (lmi, mra.project_leaf (tree_idx, element, func));
    }
  }

  for (const auto &lmi : mra.refinement_set[level])
    lmi_map->erase (lmi);

  clear_state (mra);

  return num_marked;
}

/**
 * @brief Adaptive bottom-up initialization on given initial data
 *
 * Projects onto the uniform level-1 forest, then per level thresholds the
 * details and refines the significant leaves one further level by direct
 * projection; jumping families are kept regardless. Never builds the uniform
 * max_level grid.
 */
template <typename TMultiscale, typename TFunc, typename TCriterion = hard_thresholding>
  requires coarsening_criterion<TCriterion, TMultiscale>
void
initialize_data_adaptive (TMultiscale &mra, t8_cmesh_t mesh, const t8_scheme *scheme, int max_level, TFunc &&func,
                          TCriterion criterion = {})
{
  auto c_thresh = 1.0;
  if constexpr (requires { criterion.c_thresh; })
    c_thresh = criterion.c_thresh;

  mra.initialize_data (mesh, scheme, 1, func);

  for (auto l = 1; l < max_level; ++l) {
    const auto jumps = detect_jumps (mra, l, c_thresh);

    coarsen (mra, std::max (l - 1, 1), l, jump_guarded<TCriterion, TMultiscale> { criterion, jumps });
    project_onto_children (mra, l, func);
  }
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
