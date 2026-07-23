#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/adapt/grading.hxx"
#include "t8_mra/criteria/coarsening_criterion.hxx"
#include "t8_mra/data/element_data.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_cmesh/t8_cmesh.h"
#include "t8_forest/t8_forest_general.h"
#include <t8.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace t8_mra::adapt
{

/**
 * @brief Destructive fine->coarse sweep collapsing non-significant families
 *
 * Each complete family is two-scale transformed; a non-significant one is
 * collapsed in place (parent replaces the children in lmi_map, children marked).
 * Collapsing before the sweep descends makes the parent a leaf for the next
 * level, so one traversal captures the whole cascade.
 *
 * @return Number of families collapsed on this rank
 */
template <typename MS, typename Criterion>
unsigned int
coarsen_sweep (MS &mra, int min_level, int max_level, Criterion &criterion)
{
  using element_t = typename MS::element_t;
  using detail_t = typename MS::detail_t;
  using levelmultiindex = typename MS::levelmultiindex;

  auto *lmi_map = mra.get_lmi_map ();
  auto num_marked = 0u;

  for (auto l = max_level; l > min_level; --l) {
    typename MS::index_set candidates;
    candidates.reserve (lmi_map->size (l));
    for (const auto &[lmi, _] : (*lmi_map)[l])
      candidates.insert (t8_mra::parent_lmi (lmi));

    for (const auto &parent : candidates) {
      const auto siblings = t8_mra::children_lmi (parent);

      std::array<element_t, levelmultiindex::NUM_CHILDREN> data_on_siblings;
      auto family_complete = true;

      for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
        const auto *sibling = lmi_map->find (siblings[k]);
        if (sibling == nullptr) {
          family_complete = false;
          break;
        }
        data_on_siblings[k] = *sibling;
      }

      if (!family_complete)
        continue;

      detail_t data_on_coarse;
      mra.transform.two_scale_family (data_on_siblings, data_on_coarse);
      mra.d_map.insert (parent, data_on_coarse);

      if (criterion.significant (mra, parent))
        continue;

      lmi_map->insert (parent, static_cast<const element_t &> (data_on_coarse));
      for (const auto &child : siblings) {
        lmi_map->erase (child);
        mra.coarsening_set.insert (child);
      }

      ++num_marked;
    }
  }

  return num_marked;
}

/**
 * @brief Adaptive coarsening from max_level down to min_level
 *
 * One destructive fine->coarse sweep on the maps per pass; across ranks an outer
 * fixpoint (adapt + repartition make seam families whole) until no rank marks.
 */
template <typename MS, typename Criterion = hard_thresholding>
  requires coarsening_criterion<Criterion, MS>
void
coarsen (MS &mra, int min_level, int max_level, Criterion criterion = {})
{
  if constexpr (criterion_has_prepare<Criterion, MS>)
    criterion.prepare (mra);

  for (auto pass = 0;; ++pass) {
    clear_state (mra);

    const auto num_marked = coarsen_sweep (mra, min_level, max_level, criterion);

    t8_debugf ("MRA coarsen pass %d: %u families marked, %zu leaves remain\n", pass, num_marked,
               mra.get_lmi_map ()->size ());

    if (mra.grid.global_num_marks (num_marked) == 0)
      break;

    mra.grid.adapt (MS::static_coarsening_callback, 1);
    mra.grid.repartition ();
  }

  clear_state (mra);
}

/// Per-component max mean magnitude over a level's leaves, reduced across ranks
/// (floored at 1) so every rank normalizes jump detection identically.
template <typename MS>
auto
global_v_max (MS &mra, int level)
{
  std::array<double, MS::U_DIM> local;
  local.fill (1.0);
  for (const auto &[lmi, _] : (*mra.get_lmi_map ())[level]) {
    const auto m = mra.mean_val (lmi);
    for (auto u = 0u; u < MS::U_DIM; ++u)
      local[u] = std::max (local[u], std::abs (m[u]));
  }

  std::array<double, MS::U_DIM> global;
  sc_MPI_Allreduce (local.data (), global.data (), MS::U_DIM, sc_MPI_DOUBLE, sc_MPI_MAX, mra.grid.comm);
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
template <typename MS>
auto
detect_jumps (MS &mra, int level, double c_thresh)
{
  using levelmultiindex = typename MS::levelmultiindex;

  mra.grid.ghost_exchange ();

  auto *lmi_map = mra.get_lmi_map ();
  auto &ghost_map = mra.grid.ghost_map;
  const auto v_max = global_v_max (mra, level);

  std::unordered_map<size_t, double> face_jump;
  mra.grid.for_each_face_neigh (
    [&] (const auto &lmi) { return lmi.level () == static_cast<unsigned int> (level); },
    [&] (const auto &lmi, t8_eclass_t, t8_gloidx_t, t8_element_t *, const auto &neigh_lmi) {
      const auto *neigh_data = lmi_map->contains (neigh_lmi)    ? &lmi_map->get (neigh_lmi)
                               : ghost_map.contains (neigh_lmi) ? &ghost_map.get (neigh_lmi)
                                                                : nullptr;
      if (neigh_data == nullptr)
        return;

      const auto mean_inner = mra.mean_val (lmi);
      const auto mean_neigh = mra.mean_val (*neigh_data);
      auto &diff = face_jump[lmi.index];
      for (auto u = 0u; u < MS::U_DIM; ++u)
        diff = std::max (diff, std::abs (mean_inner[u] - mean_neigh[u]) / v_max[u]);
    });

  typename MS::index_set jumps;
  for (const auto &[index, diff] : face_jump) {
    const auto lmi = levelmultiindex (index);
    const auto h = std::pow (lmi_map->get (lmi).vol, 1.0 / MS::DIM);
    if (diff > c_thresh * std::sqrt (h))
      jumps.insert (t8_mra::parent_lmi (lmi));
  }

  mra.grid.globalize (jumps);
  return jumps;
}

/// Coarsening criterion wrapper: families with a detected jump are always
/// significant. Only used by the bottom-up initialization.
template <typename Criterion, typename MS>
struct jump_guarded
{
  Criterion &criterion;
  const typename MS::index_set &jumps;

  void
  prepare (MS &mra)
  {
    if constexpr (criterion_has_prepare<Criterion, MS>)
      criterion.prepare (mra);
  }

  bool
  significant (MS &mra, const typename MS::levelmultiindex &lmi)
  {
    return jumps.contains (lmi) || criterion.significant (mra, lmi);
  }
};

/**
 * @brief Refine every leaf at the given level and project the initial data
 *
 * Unlike refine(), the children data is projected directly from the initial
 * data (exact up to quadrature), not predicted. Building block of the bottom-up
 * initialization.
 *
 * @return Number of leaves refined
 */
template <typename MS, typename Func>
unsigned int
project_onto_children (MS &mra, int level, Func &&func)
{
  clear_state (mra);

  for (const auto &[lmi, _] : (*mra.get_lmi_map ())[level])
    mra.refinement_set.insert (lmi);

  const auto num_marked = mra.refinement_set[level].size ();
  if (mra.grid.global_num_marks (num_marked) == 0)
    return 0;

  mra.grid.adapt (MS::static_refinement_callback);

  auto *lmi_map = mra.get_lmi_map ();
  auto *user_data = mra.get_user_data ();

  const auto num_local_trees = t8_forest_get_num_local_trees (mra.grid.get_forest ());
  t8_locidx_t current_idx = 0;
  for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (mra.grid.get_forest (), tree_idx);
    for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
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
template <typename MS, typename Func, typename Criterion = hard_thresholding>
  requires coarsening_criterion<Criterion, MS>
void
initialize_data_adaptive (MS &mra, t8_cmesh_t mesh, const t8_scheme *scheme, int max_level, Func &&func,
                          Criterion criterion = {})
{
  auto c_thresh = 1.0;
  if constexpr (requires { criterion.c_thresh; })
    c_thresh = criterion.c_thresh;

  mra.initialize_data (mesh, scheme, 1, func);

  for (auto l = 1; l < max_level; ++l) {
    const auto jumps = detect_jumps (mra, l, c_thresh);
    coarsen (mra, std::max (l - 1, 1), l, jump_guarded<Criterion, MS> { criterion, jumps });
    project_onto_children (mra, l, func);
  }
}

}  // namespace t8_mra::adapt

#endif  // T8_ENABLE_MRA
