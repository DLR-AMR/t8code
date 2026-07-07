#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/base.hxx"
#include "t8_mra/criteria/coarsening_criterion.hxx"
#include "t8_mra/criteria/refinement_criterion.hxx"
#include "t8_forest/t8_forest_adapt.h"
#include "t8_forest/t8_forest_ghost.h"
#include "t8_forest/t8_forest_partition.h"

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace t8_mra
{

/**
 * @brief Adaptation mixin for multiscale classes
 *
 * Provides the t8code forest adaptation functionality:
 *   - coarsen(): remove non-significant families
 *   - refine(): refine steep families and grade their neighbourhood
 *
 * Both work the same way: the full multi-level decision is computed on the
 * maps (non-destructive detail analysis, marks, data moves, iterated to
 * the cascade/grading fixpoint), then ONE recursive forest adapt realizes
 * all marks — t8code re-feeds newly formed elements to the callback within
 * the pass, so the forest is rebuilt exactly once.
 *   - the adapt machinery keeping lmi_idx and lmi_map in sync with the
 *     forest leaves
 *
 * What counts as "significant" or in need of refinement is decided by
 * exchangeable criteria: a coarsening_criterion (default: hard_thresholding,
 * see criteria/coarsening_criterion.hxx) and a refinement_criterion (default:
 * harten_prediction, see criteria/refinement_criterion.hxx).
 *
 * @tparam Derived CRTP-derived class (multiscale implementation)
 */
template <typename Derived>
class multiscale_adaptation {
 protected:
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

  /**
   * @brief Reset all per-pass multiscale state
   *
   * Stale td_set/d_map entries would corrupt later thresholding, stale
   * refinement/coarsening sets would adapt the wrong cells.
   */
  void
  clear_multiscale_state ()
  {
    derived ().d_map.erase_all ();
    derived ().td_set.erase_all ();
    derived ().refinement_set.erase_all ();
    derived ().coarsening_set.erase_all ();
  }

  /**
   * @brief Global mark count across all ranks
   *
   * t8_forest_new_adapt is collective: the decision to adapt or skip must
   * be unanimous. Adapting with an empty local set is a no-op, so ranks
   * without local marks simply participate.
   */
  unsigned int
  global_num_marks (unsigned int local_marks)
  {
    auto global_marks = local_marks;
    sc_MPI_Allreduce (&local_marks, &global_marks, 1, sc_MPI_UNSIGNED, sc_MPI_MAX, derived ().comm);
    return global_marks;
  }

  /**
   * @brief Iterate over every (leaf, same-level face neighbour) pair on this rank
   *
   * Domain-boundary faces are skipped. The visitor receives the source leaf
   * lmi, the neighbour's tree class, global tree id and element, and the
   * neighbour lmi. The neighbour element is scratch reused across a tree's
   * faces — the visitor must not retain it.
   */
  template <typename Func>
  void
  for_each_face_neigh (Func &&func)
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    auto *user_data = derived ().get_user_data ();
    auto *forest = derived ().forest;
    const auto *scheme = t8_forest_get_scheme (forest);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = 0u;

    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      t8_element_t *neigh_element;
      scheme->element_new (tree_class, 1, &neigh_element);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto num_faces = scheme->element_get_num_faces (tree_class, element);

        for (auto face = 0; face < num_faces; ++face) {
          int neigh_face;
          const auto neigh_gtreeid
            = t8_forest_element_face_neighbor (forest, tree_idx, element, neigh_element, tree_class, face, &neigh_face);

          if (neigh_gtreeid < 0)
            continue;

          func (lmi, tree_class, neigh_gtreeid, neigh_element, levelmultiindex (neigh_gtreeid, neigh_element, scheme));
        }
      }

      scheme->element_destroy (tree_class, 1, &neigh_element);
    }
  }

  /**
   * @brief Adapt the forest with the given callback and rebuild the lmi index
   *
   * Performs one t8_forest_new_adapt pass, moves the lmi_map (already
   * updated by the MST operations) into fresh user data of the new forest,
   * recomputes the per-leaf lmi index directly from the new leaves, releases
   * the old forest and runs the element-specific post-adaptation hook.
   */
  void
  adapt_forest (t8_forest_adapt_t adapt_callback, int recursive = 0)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    t8_forest_ref (derived ().forest);
    t8_forest_t new_forest
      = t8_forest_new_adapt (derived ().forest, adapt_callback, recursive, 0, derived ().get_user_data ());

    t8_mra::forest_data<element_t> *new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
    std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

    const auto num_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (new_forest);
    new_user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);
    new_user_data->mra_instance = &derived ();

    t8_forest_set_user_data (new_forest, new_user_data);

    // The LMI is computable from the element alone, so the index needs no
    // old/new leaf correspondence (cf. iterate_replace): one sweep over the
    // new leaves.
    const auto *scheme = t8_forest_get_scheme (new_forest);
    const auto num_local_trees = t8_forest_get_num_local_trees (new_forest);
    auto current_idx = t8_locidx_t { 0 };
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto gtreeid = t8_forest_global_tree_id (new_forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (new_forest, tree_idx);
      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (new_forest, tree_idx, ele_idx);
        t8_mra::set_lmi_forest_data (new_user_data, current_idx, levelmultiindex (gtreeid, element, scheme));
      }
    }

    auto *old_user_data = derived ().get_user_data ();
    delete old_user_data->lmi_map;
    sc_array_destroy (old_user_data->lmi_idx);
    T8_FREE (old_user_data);
    t8_forest_unref (&derived ().forest);

    derived ().forest = new_forest;
    // Ghost layer and ghost data do not survive the rebuild
    derived ().ghost_map.erase_all ();
    derived ().post_adaptation_hook ();
  }

 public:
  //=============================================================================
  // Adaptation Callbacks
  //=============================================================================

  /**
   * @brief Coarsening callback for t8code
   *
   * A family is coarsened iff one of its members is marked in
   * coarsening_set (all children of a pruned family are marked, so testing
   * the first member suffices regardless of child ordering). The LMI is
   * computed from the element itself, not looked up in forest_from's index:
   * this keeps the callback valid for recursive adaptation, where t8code
   * passes newly formed families that have no forest_from index.
   *
   * @return -1 to coarsen the family, 0 to keep
   */
  int
  coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, int is_family, int num_elements,
                       t8_element_t *elements[])
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    if (!is_family)
      return 0;

    const auto gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
    const auto lmi = levelmultiindex (gtreeid, elements[0], scheme);

    return derived ().coarsening_set.contains (lmi) ? -1 : 0;
  }

  /**
   * @brief Refinement callback for t8code
   *
   * A leaf is refined iff it is marked in refinement_set. The LMI is
   * computed from the element itself, not looked up in forest_from's index:
   * this keeps the callback valid for recursive adaptation, where t8code
   * passes newly created children that have no forest_from index.
   *
   * @return 1 to refine, 0 to keep
   */
  int
  refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, int is_family, int num_elements,
                       t8_element_t *elements[])
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    const auto gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
    const auto lmi = levelmultiindex (gtreeid, elements[0], scheme);

    return derived ().refinement_set.contains (lmi) ? 1 : 0;
  }

  //=============================================================================
  // Static Callback Wrappers (Required by t8code C API)
  //=============================================================================

  /**
   * @brief Static coarsening callback wrapper, routes via forest user data
   */
  static int
  static_coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;
    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    auto *mra = user_data->mra_instance;
    return static_cast<Derived *> (mra)->coarsening_callback (forest, forest_from, which_tree, tree_class,
                                                              local_ele_idx, scheme, is_family, num_elements, elements);
  }

  /**
   * @brief Static refinement callback wrapper, routes via forest user data
   */
  static int
  static_refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;
    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    auto *mra = user_data->mra_instance;
    return static_cast<Derived *> (mra)->refinement_callback (forest, forest_from, which_tree, tree_class,
                                                              local_ele_idx, scheme, is_family, num_elements, elements);
  }

  //=============================================================================
  // Coarsening
  //=============================================================================

  /**
   * @brief Perform adaptive coarsening from max_level down to min_level
   *
   * Analysis phase, a map-side cascade to the fixpoint (forest untouched):
   *   1. Compute details of all complete leaf families with parents in
   *      [min_level, max_level) via the non-destructive multiscale transformation.
   *   2. The children of a non-significant family are marked; the parent
   *      takes its data straight from the analysis and the children leave
   *      lmi_map. Significant families are never touched — no
   *      forward/inverse round-trip on kept data.
   *   A new parent can complete a coarser leaf family, which the next
   *   cascade round analyses (lmi_map IS the virtual leaf set), so
   *   multi-level coarsening resolves entirely on the maps. The cascade
   *   depth bounds the loop by max_level - min_level rounds.
   * Adaptation phase: ONE recursive forest adapt realizes all accumulated
   * marks; t8code re-feeds newly formed families to the callback within the
   * pass, collapsing the cascade into a single forest rebuild.
   *
   * @param min_level Minimum level to coarsen to
   * @param max_level Maximum level to start from
   * @param criterion Coarsening criterion (default: hard thresholding)
   */
  template <typename Criterion = hard_thresholding>
    requires coarsening_criterion<Criterion, Derived>
  void
  coarsen (int min_level, int max_level, Criterion criterion = {})
  {
    if constexpr (criterion_has_prepare<Criterion, Derived>)
      criterion.prepare (derived ());

    // Outer fixpoint: each pass coarsens every family that is complete on its
    // rank (the inner cascade peels multiple levels map-side), then adapts and
    // repartitions. The repartition is coarsening-aligned (set_for_coarsening),
    // so a family split across a rank seam — which multiscale_transformation skips as
    // incomplete — becomes whole on the next pass and finally coarsens. Repeat
    // until no rank marks anything. Serial / power-of-2 ranks have no split
    // families, so this is one real adapt plus a cheap empty pass.
    for (auto pass = 0;; ++pass) {
      clear_multiscale_state ();

      auto num_marked = 0u;
      // Re-fetch each pass: adapt_forest/repartition delete and replace lmi_map.
      auto *lmi_map = derived ().get_lmi_map ();

      for (auto round = 0; round < max_level - min_level; ++round) {
        derived ().d_map.erase_all ();

        derived ().multiscale_transformation (min_level, max_level);

        typename Derived::index_set coarsen_parents;
        for (auto L = min_level; L < max_level; ++L) {
          for (const auto &[lmi, _] : derived ().d_map[L]) {
            if (criterion.significant (derived (), lmi))
              continue;

            coarsen_parents.insert (lmi);
            for (const auto &child : t8_mra::children_lmi (lmi))
              derived ().coarsening_set.insert (child);
          }
        }

        t8_debugf ("MRA coarsen pass %d round %d: %zu families marked\n", pass, round, coarsen_parents.size ());
        if (coarsen_parents.empty ())
          break;

        num_marked += coarsen_parents.size ();

        // Data move on the maps: parent data straight from the analysis, the
        // children leave lmi_map — the next round sees the coarsened grid.
        for (const auto &lmi : coarsen_parents) {
          lmi_map->insert (lmi, derived ().d_map.get (lmi));
          for (const auto &child : t8_mra::children_lmi (lmi))
            lmi_map->erase (child);
        }
      }

      t8_debugf ("MRA coarsen pass %d: %u families marked, %zu leaves remain\n", pass, num_marked, lmi_map->size ());

      if (global_num_marks (num_marked) == 0)
        break;

      adapt_forest (static_coarsening_callback, 1);
      repartition ();
    }

    clear_multiscale_state ();
  }

  //=============================================================================
  // Refinement
  //=============================================================================

  /// Empty realized-set stand-in for callers that realize their marks per
  /// round (balance) and never descend virtually
  struct no_marks
  {
    template <typename Lmi>
    bool
    contains (const Lmi &) const
    {
      return false;
    }
  };

  /**
   * @brief Resolve a same-level neighbour lmi to its covering leaf, mark it
   *
   * The covering leaf is found by walking up the LMI hierarchy in lmi_map,
   * then descending the path back towards the neighbour through the marks
   * in realized (treated as performed refinements). A covering leaf more
   * than level_slack levels coarser than the neighbour is pulled up one
   * level (slack 0: grading, slack 1: balance).
   *
   * @param neigh_lmi Same-level neighbour of some source leaf
   * @param min_level No marks are generated below this level
   * @param level_slack Tolerated level difference
   * @param realized Marks treated as performed refinements
   * @return 1 on a new mark, 0 if nothing to do, -1 if no covering leaf is
   *         local (finer region or remote rank — caller decides)
   */
  // Lmi/RealizedSet deduced: naming Derived's nested types in the signature
  // would require the still-incomplete Derived (CRTP).
  template <typename Lmi, typename RealizedSet>
  int
  resolve_pull_up (const Lmi &neigh_lmi, int min_level, unsigned int level_slack, const RealizedSet &realized)
  {
    auto *lmi_map = derived ().get_lmi_map ();

    auto walk = neigh_lmi;
    while (walk.level () > 0 && !lmi_map->contains (walk))
      walk = t8_mra::parent_lmi (walk);

    if (!lmi_map->contains (walk))
      return -1;

    // Marks realized in earlier rounds count as performed: descend along
    // the path to the neighbour while the covering leaf is one of them.
    while (walk.level () + level_slack < neigh_lmi.level () && realized.contains (walk)) {
      auto down = neigh_lmi;
      while (down.level () > walk.level () + 1)
        down = t8_mra::parent_lmi (down);
      walk = down;
    }

    // Covering leaf too coarse -> pull it up one level
    if (walk.level () + level_slack < neigh_lmi.level () && static_cast<int> (walk.level ()) >= min_level
        && !derived ().refinement_set.contains (walk)) {
      derived ().refinement_set.insert (walk);
      return 1;
    }

    return 0;
  }

  /**
   * @brief Ship pull-up requests to their owner ranks, resolve received ones
   *
   * A same-level neighbour whose covering leaf is not local is sent to the
   * rank owning that region (payload: the raw lmi index). The owner
   * resolves the request against its own lmi_map and inserts the mark into
   * its own refinement_set — no reply needed, the mark belongs to the
   * owner. Collective: every rank must call this once per round.
   *
   * @param outgoing Per-rank request lists
   * @param min_level Passed through to resolve_pull_up
   * @param level_slack Passed through to resolve_pull_up
   * @param realized Passed through to resolve_pull_up
   * @return Number of new LOCAL marks created by received requests
   */
  template <typename RealizedSet>
  unsigned int
  exchange_pull_up_requests (const std::vector<std::vector<size_t>> &outgoing, int min_level, unsigned int level_slack,
                             const RealizedSet &realized)
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    int mpisize;
    sc_MPI_Comm_size (derived ().comm, &mpisize);

    std::vector<int> send_counts (mpisize, 0);
    for (auto rank = 0; rank < mpisize; ++rank)
      send_counts[rank] = static_cast<int> (outgoing[rank].size ());

    std::vector<int> recv_counts (mpisize, 0);
    sc_MPI_Alltoall (send_counts.data (), 1, sc_MPI_INT, recv_counts.data (), 1, sc_MPI_INT, derived ().comm);

    std::vector<std::vector<size_t>> incoming (mpisize);
    std::vector<sc_MPI_Request> requests;
    requests.reserve (2 * mpisize);

    for (auto rank = 0; rank < mpisize; ++rank) {
      if (recv_counts[rank] > 0) {
        incoming[rank].resize (recv_counts[rank]);
        requests.emplace_back ();
        sc_MPI_Irecv (incoming[rank].data (), recv_counts[rank] * sizeof (size_t), sc_MPI_BYTE, rank, 0,
                      derived ().comm, &requests.back ());
      }
      if (send_counts[rank] > 0) {
        requests.emplace_back ();
        sc_MPI_Isend (const_cast<size_t *> (outgoing[rank].data ()), send_counts[rank] * sizeof (size_t), sc_MPI_BYTE,
                      rank, 0, derived ().comm, &requests.back ());
      }
    }
    sc_MPI_Waitall (static_cast<int> (requests.size ()), requests.data (), sc_MPI_STATUSES_IGNORE);

    auto num_new_marks = 0u;
    for (const auto &batch : incoming)
      for (const auto index : batch) {
        auto neigh_lmi = levelmultiindex {};
        neigh_lmi.index = index;
        // Not resolvable here either (finer region): nothing to do
        if (resolve_pull_up (neigh_lmi, min_level, level_slack, realized) > 0)
          ++num_new_marks;
      }

    return num_new_marks;
  }

  /**
   * @brief Harten's neighbour prediction: one grading round on the maps
   *
   * A family marked by the refinement criterion (refine_neighbours -> parent
   * in td_set) requires all its same-level neighbours to carry children too
   * (Harten reference semantics). Leaf formulation used here: every leaf of
   * a marked family constructs its same-level face neighbours and pulls the
   * covering leaf up by one level if it is coarser.
   *
   * The fixpoint for level jumps larger than one is reached by repeated
   * calls WITHOUT touching the forest: marks realized in earlier rounds
   * count as performed, so the covering leaf is resolved against lmi_map and
   * then descended along the realized path towards the neighbour. Grading
   * sources are always real forest leaves — td_set only ever holds families
   * from the analysis of the actual grid (virtual children carry zero
   * details and never enter it).
   *
   * The same-level neighbour is constructed geometrically via
   * t8_forest_element_face_neighbor (no balance or ghost layer required). A
   * neighbour whose covering leaf lives on another rank is shipped to its
   * owner (exchange_pull_up_requests), which marks its own refinement_set.
   *
   * Collective: the returned mark count is global, so all ranks run the
   * same number of grading rounds.
   *
   * @param min_level No marks are generated below this level
   * @param realized Marks of EARLIER rounds (not this one), treated as
   *                 performed refinements during the descent
   * @return GLOBAL number of new marks in this round
   */
  // RealizedSet deduced (= levelindex_set): naming Derived's nested types in
  // the signature would require the still-incomplete Derived (CRTP).
  template <typename PriorRefinements>
  unsigned int
  neighbour_prediction (int min_level, const PriorRefinements &prior_refinements)
  {
    int mpirank, mpisize;
    sc_MPI_Comm_rank (derived ().comm, &mpirank);
    sc_MPI_Comm_size (derived ().comm, &mpisize);

    std::vector<std::vector<size_t>> outgoing (mpisize);
    auto num_new_marks = 0u;

    for_each_face_neigh ([&] (const auto &lmi, t8_eclass_t tree_class, t8_gloidx_t neigh_gtreeid,
                              t8_element_t *neigh_element, const auto &neigh_lmi) {
      if (lmi.level () == 0 || !derived ().td_set.contains (t8_mra::parent_lmi (lmi)))
        return;

      const auto res = refine_covering_leaf (neigh_lmi, min_level, 0u, prior_refinements);
      if (res > 0)
        ++num_new_marks;
      else if (res < 0 && mpisize > 1) {
        const auto owner = t8_forest_element_find_owner (derived ().forest, neigh_gtreeid, neigh_element, tree_class);
        if (owner != mpirank)
          outgoing[owner].push_back (neigh_lmi.index);
      }
    });

    if (mpisize > 1)
      num_new_marks += exchange_refine_requests (outgoing, min_level, 0u, prior_refinements);

    return global_num_marks (num_new_marks);
  }

  /**
   * @brief Perform adaptive refinement from min_level up to max_level
   *
   * Analysis phase (whole grid, forest untouched):
   *   1. Compute details of ALL complete leaf families via the
   *      non-destructive two-scale transform (stored at the parent levels).
   *   2. Apply the refinement criterion to each family detail (level L):
   *        - refine_neighbours: same-level neighbours must carry children
   *          -> the family grades its neighbourhood (td_set)
   *        - refine (only for L < max_level-1): the family's children
   *          (leaves at L+1) are refined one further level
   *      One pass suffices: children created by refinement carry zero
   *      details and trigger no further criterion marks.
   *   3. Grading fixpoint on the maps: each round pulls covering leaves up
   *      by one level against the marks realized in earlier rounds. A
   *      family with a marked child stops grading — once the marks are
   *      realized it is no longer a leaf family.
   * Adaptation phase: children data of all marks is reconstructed by
   * inverse two-scale with zero details (= exact polynomial subdivision);
   * multi-level mark chains descend from a real leaf through contiguous
   * marks, so the level-ascending inverse transform always finds its
   * parents. Then ONE recursive forest adapt realizes all marks; t8code
   * re-feeds new children to the callback within the pass, collapsing the
   * grading rounds into a single forest rebuild.
   *
   * Leaves at level 0 have no parent family and are never refinement
   * candidates themselves (their families are, via their parents).
   *
   * @param min_level Minimum level to start from
   * @param max_level Maximum level to refine to
   * @param criterion Refinement criterion (default: Harten's prediction)
   */
  template <typename Criterion = harten_prediction>
    requires refinement_criterion<Criterion, Derived>
  void
  refine (int min_level, int max_level, Criterion criterion = {})
  {
    if constexpr (criterion_has_prepare<Criterion, Derived>)
      criterion.prepare (derived ());

    clear_multiscale_state ();

    derived ().multiscale_transformation (0, max_level);

    auto num_families = 0u;
    for (auto L = 0; L < max_level; ++L) {
      for (const auto &[lmi, _] : derived ().d_map[L]) {
        ++num_families;

        if (criterion.refine_neighbours (derived (), lmi))
          derived ().td_set.insert (lmi);

        if (L < max_level - 1 && criterion.refine (derived (), lmi))
          for (const auto &child : t8_mra::children_lmi (lmi))
            derived ().refinement_set.insert (child);
      }
    }

    for (auto l = 0; l < min_level; ++l)
      derived ().refinement_set.erase (l);

    auto prior_refinements = derived ().refinement_set;
    prior_refinements.erase_all ();

    for (auto round = 0;; ++round) {
      const auto new_marks = neighbour_prediction (min_level, prior_refinements);
      t8_debugf ("MRA refine grading round %d: %u new marks\n", round, new_marks);
      if (new_marks == 0)
        break;

      prior_refinements = derived ().refinement_set;

      typename Derived::index_set stopped;
      for (auto L = 0; L < max_level; ++L)
        for (const auto &lmi : derived ().td_set[L]) {
          const auto children = t8_mra::children_lmi (lmi);
          if (std::any_of (children.begin (), children.end (),
                           [&] (const auto &child) { return derived ().refinement_set.contains (child); }))
            stopped.insert (lmi);
        }
      for (const auto &lmi : stopped)
        derived ().td_set.erase (lmi);
    }

    const auto num_marked = num_refinement_marks (min_level, max_level);
    t8_debugf ("MRA refine analysis: %u leaf families, %u leaves marked\n", num_families, num_marked);

    if (global_num_marks (num_marked) == 0) {
      clear_multiscale_state ();
      return;
    }

    //--------------------------------------------------------------------------
    // Adaptation phase
    //--------------------------------------------------------------------------

    // Reconstruct children data for all marks: inverse two-scale with zero
    // details. This moves the data one level down in lmi_map per mark,
    // exactly mirroring what the forest adapt below does.
    derived ().d_map.erase_all ();
    for (auto l = min_level; l < max_level; ++l)
      for (const auto &lmi : derived ().refinement_set[l])
        derived ().d_map.insert (lmi, typename Derived::detail_t {});

    derived ().inverse_multiscale_transformation (min_level, max_level);

    // ONE recursive forest adapt realizes all marks across all levels
    adapt_forest (static_refinement_callback, 1);
    repartition ();

    clear_multiscale_state ();
  }

  //=============================================================================
  // Balancing
  //=============================================================================

  /**
   * @brief One balancing round: pull up covering leaves across faces
   *
   * Every leaf constructs its same-level face neighbours geometrically and
   * resolves them to the covering leaf via the LMI hierarchy. A covering
   * leaf more than one level coarser is marked and refined one level;
   * children data comes from the inverse two-scale transform with zero
   * details, so the represented data is unchanged.
   *
   * @return Number of leaves marked in this round
   */
  unsigned int
  balance_round ()
  {
    clear_multiscale_state ();

    int mpirank, mpisize;
    sc_MPI_Comm_rank (derived ().comm, &mpirank);
    sc_MPI_Comm_size (derived ().comm, &mpisize);
    std::vector<std::vector<size_t>> outgoing (mpisize);

    for_each_face_neigh ([&] (const auto &lmi, t8_eclass_t tree_class, t8_gloidx_t neigh_gtreeid,
                              t8_element_t *neigh_element, const auto &neigh_lmi) {
      if (lmi.level () < 2)
        return;

      if (refine_covering_leaf (neigh_lmi, 0, 1u, no_prior_marks {}) < 0 && mpisize > 1) {
        const auto owner = t8_forest_element_find_owner (derived ().forest, neigh_gtreeid, neigh_element, tree_class);
        if (owner != mpirank)
          outgoing[owner].push_back (neigh_lmi.index);
      }
    });

    if (mpisize > 1)
      exchange_refine_requests (outgoing, 0, 1u, no_prior_marks {});

    const auto num_marked = global_num_marks (num_refinement_marks (0, derived ().maximum_level));
    if (num_marked == 0) {
      clear_multiscale_state ();
      return 0;
    }

    apply_refinement (0, derived ().maximum_level, 0);
    clear_multiscale_state ();

    return num_marked;
  }

  /**
   * @brief Restore the 2:1 face balance of the grid
   *
   * Rounds iterate until no leaf has a face neighbour more than one level
   * coarser. Each round lifts violating covering leaves by one level, so a
   * jump of k levels resolves in k-1 rounds; refining a leaf can create new
   * violations against its own coarser neighbours, which the next round
   * catches. Terminates: every round refines at least one leaf and levels
   * are bounded by max_level.
   */
  void
  balance ()
  {
    auto rounds = 0;
    while (balance_round () > 0)
      t8_debugf ("MRA balance round %d\n", rounds++);

    if (rounds > 0)
      repartition ();
  }

  //=============================================================================
  // Repartitioning
  //=============================================================================

  /**
   * @brief Repartition the forest along the SFC and migrate the leaf data
   *
   * Adaptive passes unbalance the per-rank load. Builds the partitioned
   * forest, ships each leaf's element_data along the partition
   * (t8_forest_partition_data, SFC leaf order), then rebuilds lmi_map and
   * lmi_idx from the new local leaves. Requires the standing invariant
   * (forest leaves <-> lmi_map keys 1:1); collective.
   *
   * Partitions with set_for_coarsening = 1 (t8code PFC): boundaries are
   * rounded so no family of same-level siblings is split across ranks, so
   * the next coarsen pass can merge a family that was a seam family before.
   */
  void
  repartition ()
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    static_assert (std::is_trivially_copyable_v<element_t>, "element data is shipped as raw bytes");

    int mpisize;
    sc_MPI_Comm_size (derived ().comm, &mpisize);
    if (mpisize == 1)
      return;

    auto *forest = derived ().forest;
    auto *user_data = derived ().get_user_data ();
    auto *lmi_map = derived ().get_lmi_map ();

    // Pack the leaf data in SFC leaf order
    const auto num_old = t8_forest_get_local_num_leaf_elements (forest);
    auto *data_in = sc_array_new_count (sizeof (element_t), num_old);
    for (t8_locidx_t i = 0; i < num_old; ++i) {
      const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, i);
      *reinterpret_cast<element_t *> (sc_array_index (data_in, i)) = lmi_map->get (lmi);
    }

    t8_forest_ref (forest);
    t8_forest_t new_forest;
    t8_forest_init (&new_forest);
    // set_for_coarsening = 1: PFC rounds partition boundaries so no family of
    // same-level siblings is split across ranks, letting seam families coarsen.
    t8_forest_set_partition (new_forest, forest, 1);
    t8_forest_commit (new_forest);

    const auto num_new = t8_forest_get_local_num_leaf_elements (new_forest);
    auto *data_out = sc_array_new_count (sizeof (element_t), num_new);
    t8_forest_partition_data (forest, new_forest, data_in, data_out);
    sc_array_destroy (data_in);

    // Fresh user data; lmi_map and lmi_idx rebuilt from the new leaves (the
    // migrated data arrives in the new SFC leaf order)
    auto *new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
    new_user_data->lmi_idx
      = sc_array_new_count (sizeof (levelmultiindex), num_new + t8_forest_get_num_ghosts (new_forest));
    new_user_data->mra_instance = &derived ();
    t8_forest_set_user_data (new_forest, new_user_data);

    const auto *scheme = t8_forest_get_scheme (new_forest);
    const auto num_local_trees = t8_forest_get_num_local_trees (new_forest);
    auto current_idx = t8_locidx_t { 0 };
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto gtreeid = t8_forest_global_tree_id (new_forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (new_forest, tree_idx);
      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto *element = t8_forest_get_leaf_element_in_tree (new_forest, tree_idx, ele_idx);
        const auto lmi = levelmultiindex (gtreeid, element, scheme);
        t8_mra::set_lmi_forest_data (new_user_data, current_idx, lmi);
        new_user_data->lmi_map->insert (lmi, *reinterpret_cast<element_t *> (sc_array_index (data_out, current_idx)));
      }
    }
    sc_array_destroy (data_out);

    delete user_data->lmi_map;
    sc_array_destroy (user_data->lmi_idx);
    T8_FREE (user_data);
    t8_forest_unref (&derived ().forest);

    derived ().forest = new_forest;
    // Ghost layer and ghost data do not survive the rebuild
    derived ().ghost_map.erase_all ();
    derived ().post_adaptation_hook ();
  }

  //=============================================================================
  // Ghost Exchange
  //=============================================================================

  /**
   * @brief Fill the ghost layer: lmi index and ghost element data
   *
   * Creates the face-ghost layer on demand (t8_forest_ghost_create handles
   * unbalanced forests), extends lmi_idx by the ghost slots and exchanges
   * them, then ships the element data of all remote leaves touching the
   * local partition into ghost_map.
   *
   * The ghost data is a read-only snapshot for solver stencils: every
   * coarsen/refine/balance/repartition invalidates it (ghost_map is
   * cleared) — call again after adapting. Adaptation itself never needs
   * ghosts; remote grading is handled by exchange_pull_up_requests.
   * Collective.
   */
  void
  ghost_exchange ()
  {
    using element_t = typename Derived::element_t;

    auto *forest = derived ().forest;

    if (t8_forest_get_num_ghosts (forest) == 0)
      t8_forest_ghost_create (forest);

    const auto num_local = t8_forest_get_local_num_leaf_elements (forest);
    const auto num_ghosts = t8_forest_get_num_ghosts (forest);

    // lmi_idx was allocated before the ghost layer existed
    auto *user_data = derived ().get_user_data ();
    sc_array_resize (user_data->lmi_idx, num_local + num_ghosts);
    t8_forest_ghost_exchange_data (forest, user_data->lmi_idx);

    // Leaf data in leaf order; the exchange fills the ghost slots
    auto *data = sc_array_new_count (sizeof (element_t), num_local + num_ghosts);
    auto *lmi_map = derived ().get_lmi_map ();
    for (t8_locidx_t i = 0; i < num_local; ++i)
      *reinterpret_cast<element_t *> (sc_array_index (data, i))
        = lmi_map->get (t8_mra::get_lmi_from_forest_data (user_data, i));

    t8_forest_ghost_exchange_data (forest, data);

    derived ().ghost_map.erase_all ();
    for (auto i = num_local; i < num_local + num_ghosts; ++i)
      derived ().ghost_map.insert (t8_mra::get_lmi_from_forest_data (user_data, i),
                                   *reinterpret_cast<element_t *> (sc_array_index (data, i)));
    sc_array_destroy (data);
  }

  //=============================================================================
  // Bottom-up Initialization
  //=============================================================================

  /**
   * @brief Mean-value jump detection on the leaves of one level
   *
   * Compares each leaf's component means against its same-level face
   * neighbours, normalized per component by the mean magnitude where it
   * exceeds 1. Marks the leaf's family when the difference exceeds
   * c_thresh * sqrt(h): smooth data decays as O(h) and falls below, a
   * discontinuity stays O(1) and keeps firing on every level. Catches jumps
   * aligned with family boundaries, where the detail coefficients vanish.
   *
   * @param level Leaves of this level are checked
   * @param c_thresh Threshold constant
   * @return Parent lmis of the jumping families
   */
  auto
  detect_jumps (int level, double c_thresh)
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    typename Derived::index_set jumps;

    auto *user_data = derived ().get_user_data ();
    auto *forest = derived ().forest;
    const auto *scheme = t8_forest_get_scheme (forest);
    auto *lmi_map = derived ().get_lmi_map ();

    std::array<double, Derived::U_DIM> v_max;
    v_max.fill (1.0);
    for (const auto &[lmi, _] : (*lmi_map)[level]) {
      const auto m = derived ().mean_val (lmi);
      for (auto u = 0u; u < Derived::U_DIM; ++u)
        v_max[u] = std::max (v_max[u], std::abs (m[u]));
    }

    // All ranks must normalize identically
    auto v_max_local = v_max;
    sc_MPI_Allreduce (v_max_local.data (), v_max.data (), Derived::U_DIM, sc_MPI_DOUBLE, sc_MPI_MAX, derived ().comm);

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = t8_locidx_t { 0 };
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      t8_element_t *neigh_element;
      scheme->element_new (tree_class, 1, &neigh_element);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        if (lmi.level () != static_cast<unsigned int> (level))
          continue;

        const auto mean_inner = derived ().mean_val (lmi);
        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto num_faces = scheme->element_get_num_faces (tree_class, element);

        auto max_diff = 0.0;
        for (auto face = 0; face < num_faces; ++face) {
          int neigh_face;
          const auto neigh_gtreeid
            = t8_forest_element_face_neighbor (forest, tree_idx, element, neigh_element, tree_class, face, &neigh_face);

          // Domain boundary
          if (neigh_gtreeid < 0)
            continue;

          // Same-level neighbour only; coarser neighbours are skipped
          const auto neigh_lmi = levelmultiindex (neigh_gtreeid, neigh_element, scheme);
          if (!lmi_map->contains (neigh_lmi))
            continue;

          const auto mean_neigh = derived ().mean_val (neigh_lmi);
          for (auto u = 0u; u < Derived::U_DIM; ++u)
            max_diff = std::max (max_diff, std::abs (mean_inner[u] - mean_neigh[u]) / v_max[u]);
        }

        const auto h = std::pow (lmi_map->get (lmi).vol, 1.0 / Derived::DIM);
        if (max_diff > c_thresh * std::sqrt (h))
          jumps.insert (t8_mra::parent_lmi (lmi));
      }

      scheme->element_destroy (tree_class, 1, &neigh_element);
    }

    return jumps;
  }

  /**
   * @brief Coarsening criterion wrapper: families with a detected jump are
   * always significant. Only used by the bottom-up initialization.
   */
  template <typename Criterion>
  struct jump_guarded
  {
    Criterion &criterion;
    const typename Derived::index_set &jumps;

    void
    prepare (Derived &mra)
    {
      if constexpr (criterion_has_prepare<Criterion, Derived>)
        criterion.prepare (mra);
    }

    bool
    significant (Derived &mra, const typename Derived::levelmultiindex &lmi)
    {
      return jumps.contains (lmi) || criterion.significant (mra, lmi);
    }
  };

  /**
   * @brief Refine every leaf at the given level and project the initial data
   *
   * Unlike refine(), the children data is not predicted by the inverse
   * two-scale transform but projected directly from the initial data, which
   * is exact up to quadrature. Building block of the bottom-up
   * initialization.
   *
   * @param level Level whose leaves are refined (children appear at level+1)
   * @param func Initial data to project onto the new leaves
   * @return Number of leaves refined
   */
  template <typename Func>
  unsigned int
  refine_by_projection (int level, Func &&func)
  {
    clear_multiscale_state ();

    for (const auto &[lmi, _] : (*derived ().get_lmi_map ())[level])
      derived ().refinement_set.insert (lmi);

    const auto num_marked = derived ().refinement_set[level].size ();
    if (global_num_marks (num_marked) == 0)
      return 0;

    adapt_forest (static_refinement_callback);

    // The forest already carries the children, lmi_map still the parents.
    // Project the initial data onto every leaf without map entry (exactly
    // the new children), then drop the refined parents.
    auto *lmi_map = derived ().get_lmi_map ();
    auto *user_data = derived ().get_user_data ();

    const auto num_local_trees = t8_forest_get_num_local_trees (derived ().forest);
    auto current_idx = t8_locidx_t { 0 };
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (derived ().forest, tree_idx);
      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        if (lmi_map->contains (lmi))
          continue;

        const auto *element = t8_forest_get_leaf_element_in_tree (derived ().forest, tree_idx, ele_idx);
        lmi_map->insert (lmi, derived ().project_leaf (tree_idx, element, func));
      }
    }

    for (const auto &lmi : derived ().refinement_set[level])
      lmi_map->erase (lmi);

    clear_multiscale_state ();

    return num_marked;
  }

  /**
   * @brief Adaptive bottom-up initialization on given initial data
   *
   * Projects onto the uniform level-1 forest, then per level thresholds the
   * details with the coarsening criterion and refines the significant leaves
   * one further level by direct projection. Families with a mean-value jump
   * across faces are kept regardless of their details (detect_jumps). Never
   * builds the uniform max_level grid.
   *
   * @param mesh Coarse mesh
   * @param scheme Element scheme
   * @param max_level Maximum refinement level
   * @param func Initial data to project
   * @param criterion Coarsening criterion (default: hard thresholding)
   */
  template <typename Func, typename Criterion = hard_thresholding>
    requires coarsening_criterion<Criterion, Derived>
  void
  initialize_data_adaptive (t8_cmesh_t mesh, const t8_scheme *scheme, int max_level, Func &&func,
                            Criterion criterion = {})
  {
    // Jump tolerance follows the criterion's threshold constant where it has one
    auto c_thresh = 1.0;
    if constexpr (requires { criterion.c_thresh; })
      c_thresh = criterion.c_thresh;

    derived ().initialize_data (mesh, scheme, 1, func);

    for (auto l = 1; l < max_level; ++l) {
      const auto jumps = detect_jumps (l, c_thresh);
      coarsen (std::max (l - 1, 1), l, jump_guarded<Criterion> { criterion, jumps });
      refine_by_projection (l, func);
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
