#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_base.hpp"
#include "t8_mra/t8_mra_coarsening_criterion.hpp"
#include "t8_mra/t8_mra_refinement_criterion.hpp"
#include "t8_forest/t8_forest_adapt.h"
#include "t8_forest/t8_forest_iterate.h"

#include <algorithm>
#include <cmath>

namespace t8_mra
{

/**
 * @brief Adaptation mixin for multiscale classes
 *
 * Provides the t8code forest adaptation functionality:
 *   - coarsen(): remove non-significant detail information level by level
 *   - refine(): refine steep families and grade their neighbourhood
 *   - the adapt/iterate_replace machinery keeping lmi_idx and lmi_map
 *     in sync with the forest leaves
 *
 * What counts as "significant" or in need of refinement is decided by
 * exchangeable criteria: a coarsening_criterion (default: hard_thresholding,
 * see t8_mra_coarsening_criterion.hpp) and a refinement_criterion (default:
 * harten_prediction, see t8_mra_refinement_criterion.hpp).
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
   * @brief Adapt the forest with the given callback and rebuild the lmi index
   *
   * Performs one t8_forest_new_adapt pass, moves the lmi_map (already
   * updated by the MST operations) into fresh user data of the new forest,
   * rebuilds the per-leaf lmi index via iterate_replace, releases the old
   * forest and runs the element-specific post-adaptation hook.
   */
  void
  adapt_forest (t8_forest_adapt_t adapt_callback)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    t8_forest_ref (derived ().forest);
    t8_forest_t new_forest = t8_forest_new_adapt (derived ().forest, adapt_callback, 0, 0, derived ().get_user_data ());

    t8_mra::forest_data<element_t> *new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
    new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
    std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

    const auto num_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
    const auto num_ghost_elements = t8_forest_get_num_ghosts (new_forest);
    new_user_data->lmi_idx = sc_array_new_count (sizeof (levelmultiindex), num_local_elements + num_ghost_elements);
    new_user_data->mra_instance = &derived ();

    t8_forest_set_user_data (new_forest, new_user_data);
    t8_forest_iterate_replace (new_forest, derived ().forest, static_iterate_replace_callback);

    auto *old_user_data = derived ().get_user_data ();
    delete old_user_data->lmi_map;
    sc_array_destroy (old_user_data->lmi_idx);
    T8_FREE (old_user_data);
    t8_forest_unref (&derived ().forest);

    derived ().forest = new_forest;
    derived ().post_adaptation_hook ();
  }

 public:
  //=============================================================================
  // Adaptation Callbacks
  //=============================================================================

  /**
   * @brief Coarsening callback for t8code
   *
   * A family is coarsened iff its (first) member is marked in
   * coarsening_set. The stored LMI encodes the level, so set membership
   * is exact.
   *
   * @return -1 to coarsen the family, 0 to keep
   */
  int
  coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, int is_family, int num_elements,
                       t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;

    if (!is_family)
      return 0;

    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    const auto offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_ele_idx + offset);

    return derived ().coarsening_set.contains (lmi) ? -1 : 0;
  }

  /**
   * @brief Refinement callback for t8code
   *
   * A leaf is refined iff it is marked in refinement_set. The stored LMI
   * encodes the level, so set membership is exact.
   *
   * @return 1 to refine, 0 to keep
   */
  int
  refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, int is_family, int num_elements,
                       t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;

    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    const auto offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_ele_idx + offset);

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
  // Replace Callback (rebuilds the per-leaf lmi index)
  //=============================================================================

  /**
   * @brief Replace callback: updates the lmi_idx array to match the new forest
   *
   * Only the lmi_idx array is updated here. The lmi_map was already
   * populated by the MST operations and swapped to the new user data.
   */
  void
  iterate_replace_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                            const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                            t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    auto *old_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_old));
    auto *new_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_new));

    // Adjust indices with tree offsets
    const auto first_incoming_in_tree = first_incoming;
    first_incoming += t8_forest_get_tree_element_offset (forest_new, which_tree);
    first_outgoing += t8_forest_get_tree_element_offset (forest_old, which_tree);

    const auto old_lmi = t8_mra::get_lmi_from_forest_data (old_user_data, first_outgoing);

    if (refine == 0) {
      // No change: copy LMI as-is
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, old_lmi);
    }
    else if (refine == -1) {
      // Coarsening: set parent LMI
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, t8_mra::parent_lmi (old_lmi));
    }
    else {
      // Refinement: construct each child's LMI from the actual incoming
      // element. The forest's child ordering (Bey, type-dependent for
      // triangles) does not generally match the LMI reference ordering, so
      // positional assignment of children_lmi(old_lmi) would mislabel cells.
      const auto gtreeid = t8_forest_global_tree_id (forest_new, which_tree);
      for (int i = 0; i < num_incoming; ++i) {
        const auto *child_element
          = t8_forest_get_leaf_element_in_tree (forest_new, which_tree, first_incoming_in_tree + i);
        const auto child_lmi = levelmultiindex (gtreeid, child_element, scheme);
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, child_lmi);
      }
    }
  }

  /**
   * @brief Static wrapper for iterate_replace_callback, routes via forest user data
   */
  static void
  static_iterate_replace_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                   const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                                   t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    using element_t = typename Derived::element_t;
    auto *old_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_old));

    auto *mra = old_user_data->mra_instance;
    static_cast<Derived *> (mra)->iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme,
                                                            refine, num_outgoing, first_outgoing, num_incoming,
                                                            first_incoming);
  }

  //=============================================================================
  // Coarsening
  //=============================================================================

  /**
   * @brief Perform adaptive coarsening from max_level down to min_level
   *
   * Per level l (fine to coarse, levels depend on each other):
   *   1. Two-scale transform of all complete leaf families at level l
   *      (parent data + details at level l-1).
   *   2. Apply the criterion; non-significant families lose their details
   *      and their children (the current leaves) are marked for coarsening.
   *   3. Inverse transform reconstructs the children of significant families
   *      exactly (details were kept).
   *   4. Forest adapt removes the marked families.
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

    for (auto l = max_level; l > min_level; --l) {
      clear_multiscale_state ();

      derived ().multiscale_transformation (l - 1, l);

      // Prune non-significant families: drop their details and mark their
      // children for coarsening. Iterate over a copy since we erase entries.
      const auto details = derived ().d_map[l - 1];
      for (const auto &[lmi, _] : details) {
        if (!criterion.significant (derived (), lmi)) {
          derived ().d_map.erase (lmi);
          for (const auto &child : t8_mra::children_lmi (lmi))
            derived ().coarsening_set.insert (child);
        }
      }

      derived ().inverse_multiscale_transformation (l - 1, l);

      t8_debugf ("MRA coarsen pass %d: coarsening %zu of %zu leaves\n", l, derived ().coarsening_set[l].size (),
                 derived ().get_user_data ()->lmi_map->size ());

      adapt_forest (static_coarsening_callback);
    }

    clear_multiscale_state ();
  }

  //=============================================================================
  // Refinement
  //=============================================================================

  /**
   * @brief Harten's neighbour prediction on the leaf level
   *
   * A family marked by the refinement criterion (refine_neighbours -> parent
   * in td_set) requires all its same-level neighbours to carry children too
   * (Harten reference semantics). Leaf formulation used here: every leaf of
   * a marked family constructs its same-level face neighbours and pulls the
   * covering leaf up by one level if it is coarser.
   * Repeated refine calls reach the fixpoint for level jumps larger
   * than one.
   *
   * The same-level neighbour is constructed geometrically via
   * t8_forest_element_face_neighbor (no balance or ghost layer required) and
   * resolved to the covering leaf by walking up the LMI hierarchy in lmi_map.
   *
   * TODO MPI: a neighbour lmi whose covering leaf is not local resolves to
   * "not found" here. Collect those lmis per owner rank and exchange them
   * (cf. reference implementation: data_to_send), the owner then resolves
   * them against its local lmi_map.
   */
  void
  neighbour_prediction ()
  {
    using levelmultiindex = typename Derived::levelmultiindex;

    auto *user_data = derived ().get_user_data ();
    auto *forest = derived ().forest;
    const auto *scheme = t8_forest_get_scheme (forest);
    auto *lmi_map = derived ().get_lmi_map ();

    const auto num_local_trees = t8_forest_get_num_local_trees (forest);
    auto current_idx = t8_locidx_t { 0 };
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      // Scratch element for the constructed neighbours, one per tree
      t8_element_t *neigh_element;
      scheme->element_new (tree_class, 1, &neigh_element);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        if (lmi.level () == 0)
          continue;

        // Only leaves of marked families (refine_neighbours) grade their neighbourhood
        if (!derived ().td_set.contains (t8_mra::parent_lmi (lmi)))
          continue;

        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto num_faces = scheme->element_get_num_faces (tree_class, element);

        for (auto face = 0; face < num_faces; ++face) {
          int neigh_face;
          const auto neigh_gtreeid
            = t8_forest_element_face_neighbor (forest, tree_idx, element, neigh_element, tree_class, face, &neigh_face);

          // Domain boundary
          if (neigh_gtreeid < 0)
            continue;

          // Same-level neighbour as LMI, then resolve to the covering leaf:
          // walk up until an entry exists in lmi_map.
          auto walk = levelmultiindex (neigh_gtreeid, neigh_element, scheme);
          while (walk.level () > 0 && !lmi_map->contains (walk))
            walk = t8_mra::parent_lmi (walk);

          // Covering leaf coarser than this leaf -> pull it up one level.
          // Not found at all: neighbour region is refined finer (children
          // handle themselves) or lives on another rank (TODO MPI above).
          if (lmi_map->contains (walk) && walk.level () < lmi.level ())
            derived ().refinement_set.insert (walk);
        }
      }

      scheme->element_destroy (tree_class, 1, &neigh_element);
    }
  }

  /**
   * @brief One refinement round: analysis, then a single forest adapt
   *
   * Analysis phase (whole grid, before any adaptation):
   *   1. Compute details of ALL complete leaf families via the
   *      non-destructive two-scale transform (stored at the parent levels).
   *   2. Apply the refinement criterion to each family detail (level L):
   *        - refine_neighbours: same-level neighbours must carry children
   *          -> neighbour_prediction pulls coarser neighbour leaves up
   *        - refine (only for L < max_level-1): the family's children
   *          (leaves at L+1) are refined one further level
   * Adaptation phase: children data of all marked leaves is reconstructed by
   * inverse two-scale with zero details (= exact polynomial subdivision),
   * then a single forest adapt refines exactly the marked leaves; lmi_map
   * and forest stay in sync by construction.
   *
   * Leaves at level 0 have no parent family and are never refinement
   * candidates themselves (their families are, via their parents).
   *
   * @param min_level Minimum level to start from
   * @param max_level Maximum level to refine to
   * @param criterion Refinement criterion
   * @return Number of leaves marked for refinement in this round
   */
  template <typename Criterion>
    requires refinement_criterion<Criterion, Derived>
  unsigned int
  refine_round (int min_level, int max_level, Criterion &criterion)
  {
    using element_t = typename Derived::element_t;

    clear_multiscale_state ();

    //--------------------------------------------------------------------------
    // Analysis phase
    //--------------------------------------------------------------------------

    // 1. Details of all complete leaf families; lmi_map stays untouched.
    derived ().compute_leaf_details (0, max_level);

    // 2. Apply the criterion to every family detail
    auto num_families = 0u;

    for (auto L = 0; L < max_level; ++L) {
      for (const auto &[lmi, _] : derived ().d_map[L]) {
        ++num_families;

        // Remember families whose neighbourhood must be graded
        if (criterion.refine_neighbours (derived (), lmi))
          derived ().td_set.insert (lmi);

        // Refine the family's children (leaves at L+1).
        // Guard keeps the result within max_level.
        if (L < max_level - 1 && criterion.refine (derived (), lmi))
          for (const auto &child : t8_mra::children_lmi (lmi))
            derived ().refinement_set.insert (child);
      }
    }

    // 3. Grade the neighbourhood of the marked families
    neighbour_prediction ();

    // Details are consumed; only the refinement_set is needed from here on
    derived ().d_map.erase_all ();

    // Marks below min_level are outside the requested range
    for (auto l = 0; l < min_level; ++l)
      derived ().refinement_set.erase (l);

    auto num_marked = 0u;
    for (auto l = min_level; l < max_level; ++l)
      num_marked += derived ().refinement_set[l].size ();
    t8_debugf ("MRA refine analysis: %u leaf families, %zu grading neighbourhoods, %u leaves marked\n", num_families,
               derived ().td_set.size (), num_marked);

    if (num_marked == 0) {
      clear_multiscale_state ();
      return 0;
    }

    //--------------------------------------------------------------------------
    // Adaptation phase
    //--------------------------------------------------------------------------

    // Reconstruct children data for all marked leaves: inverse two-scale
    // with zero details. This moves each refined leaf's data one level down
    // in lmi_map, exactly mirroring what the forest adapt below does.
    for (auto l = min_level; l < max_level; ++l)
      for (const auto &lmi : derived ().refinement_set[l])
        derived ().d_map.insert (lmi, element_t {});

    derived ().inverse_multiscale_transformation (min_level, max_level);

    // A single forest adapt refines all marked leaves across all levels
    adapt_forest (static_refinement_callback);

    clear_multiscale_state ();

    return num_marked;
  }

  /**
   * @brief Perform adaptive refinement from min_level up to max_level
   *
   * Iterates refinement rounds until the grading fixpoint is reached:
   * the neighbour prediction lifts a coarser neighbour of a marked family
   * by one level per round, so a neighbour across a multi-level jump needs
   * several rounds to arrive at the family's leaf level. Newly created
   * children carry zero details and trigger no further marks, hence the
   * iteration terminates (bounded by max_level rounds).
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

    for (auto round = 0; round < max_level; ++round) {
      t8_debugf ("MRA refine round %d\n", round);
      if (refine_round (min_level, max_level, criterion) == 0)
        break;
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
