#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_base.hpp"
#include "t8_forest/t8_forest_adapt.h"
#include "t8_forest/t8_forest_ghost.h"
#include "t8_forest/t8_forest_iterate.h"

#include <algorithm>
#include <cmath>

namespace t8_mra
{

/**
 * @brief Adaptation mixin for multiscale classes
 *
 * Provides t8code forest adaptation functionality:
 *   - Coarsening callbacks
 *   - Refinement callbacks
 *   - Forest adaptation loops
 *   - Ghost exchange
 *   - Balancing
 *
 * @tparam Derived CRTP-derived class (multiscale implementation)
 */
template <typename Derived>
class multiscale_adaptation {
 protected:
  //=============================================================================
  // Helper: Access derived class
  //=============================================================================

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

 public:
  //=============================================================================
  // Adaptation Callbacks
  //=============================================================================

  /**
   * @brief Coarsening callback for t8code
   *
   * Checks if element is marked for coarsening in coarsening_set.
   *
   * @return -1 to coarsen family, 0 to keep
   */
  int
  coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                       t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, int is_family, int num_elements,
                       t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    if (!is_family)
      return 0;

    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    // Only process elements at current refinement level
    const auto element_level = scheme->element_get_level (tree_class, elements[0]);
    if (element_level != user_data->current_refinement_level)
      return 0;

    // Get LMI for the first child from forest_data (same approach as old version)
    const auto offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    const auto elem_idx = local_ele_idx + offset;

    const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, elem_idx);

    return derived ().coarsening_set.contains (lmi) ? -1 : 0;
  }

  /**
   * @brief Refinement callback for t8code
   *
   * Checks if element is marked for refinement in refinement_set.
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

    // Only process elements at current refinement level
    const auto element_level = scheme->element_get_level (tree_class, elements[0]);
    if (element_level != user_data->current_refinement_level)
      return 0;

    // Read the stored LMI (same approach as coarsening_callback) instead of
    // reconstructing it from the element
    const auto offset = t8_forest_get_tree_element_offset (forest_from, which_tree);
    const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_ele_idx + offset);

    return derived ().refinement_set.contains (lmi) ? 1 : 0;
  }

  //=============================================================================
  // Static Callback Wrappers (Required by t8code C API)
  //=============================================================================

  /**
   * @brief Static coarsening callback wrapper
   */
  static int
  static_coarsening_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;
    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    // Get the multiscale object from user data
    auto *mra = user_data->mra_instance;
    return static_cast<Derived *> (mra)->coarsening_callback (forest, forest_from, which_tree, tree_class,
                                                              local_ele_idx, scheme, is_family, num_elements, elements);
  }

  /**
   * @brief Static refinement callback wrapper
   */
  static int
  static_refinement_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_eclass_t tree_class, t8_locidx_t local_ele_idx, const t8_scheme_c *scheme,
                              int is_family, int num_elements, t8_element_t *elements[])
  {
    using element_t = typename Derived::element_t;
    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    // Get the multiscale object from user data
    auto *mra = user_data->mra_instance;
    return static_cast<Derived *> (mra)->refinement_callback (forest, forest_from, which_tree, tree_class,
                                                              local_ele_idx, scheme, is_family, num_elements, elements);
  }

  //=============================================================================
  // Replace Callback (For copying data between old and new forest)
  //=============================================================================

  /**
   * @brief Replace callback: updates lmi_idx array to match forest structure
   *
   * IMPORTANT: This callback should ONLY update the lmi_idx array (via set_lmi_forest_data).
   * The lmi_map was already populated correctly by the MST operations and swapped to new_user_data.
   * DO NOT modify lmi_map here!
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
      const auto parent_lmi = t8_mra::parent_lmi (old_lmi);
      t8_mra::set_lmi_forest_data (new_user_data, first_incoming, parent_lmi);
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
   * @brief Static wrapper for iterate_replace_callback
   */
  static void
  static_iterate_replace_callback (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                   const t8_eclass_t tree_class, const t8_scheme *scheme, int refine, int num_outgoing,
                                   t8_locidx_t first_outgoing, int num_incoming, t8_locidx_t first_incoming)
  {
    using element_t = typename Derived::element_t;
    auto *old_user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_old));

    // Get the multiscale object from user data
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
   * @param min_level Minimum level to coarsen to
   * @param max_level Maximum level to start from
   */
  void
  coarsen (int min_level, int max_level)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    /// Scaling due to (2.39)
    derived ().c_scaling = derived ().threshold_scaling_factor ();

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (derived ().forest);

      derived ().get_user_data ()->current_refinement_level = l;

      // Start each pass with clean multiscale state; stale td_set/d_map
      // entries would corrupt the thresholding (sync_d_with_td inserts
      // empty elements for td entries missing from d_map)
      derived ().d_map.erase_all ();
      derived ().td_set.erase_all ();
      derived ().refinement_set.erase_all ();
      derived ().coarsening_set.erase_all ();

      derived ().multiscale_transformation (l - 1, l);
      derived ().threshold (l - 1, l);
      // restore_balancing (l - 1, l); /// TODO
      derived ().generate_td_tree (l - 1, l);
      derived ().sync_d_with_td (min_level, max_level);

      derived ().inverse_multiscale_transformation (l - 1, l);

      t8_debugf ("MRA coarsen pass %d: coarsening %zu of %zu leaves\n", l, derived ().coarsening_set[l].size (),
                 derived ().get_user_data ()->lmi_map->size ());

      new_forest
        = t8_forest_new_adapt (derived ().forest, static_coarsening_callback, 0, 0, derived ().get_user_data ());

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
      std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      new_user_data->mra_instance = &derived ();
      new_user_data->current_refinement_level = l;

      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (new_forest, derived ().forest, static_iterate_replace_callback);

      // Cleanup old forest and user data
      auto *old_user_data = derived ().get_user_data ();
      delete old_user_data->lmi_map;
      sc_array_destroy (old_user_data->lmi_idx);
      T8_FREE (old_user_data);
      t8_forest_unref (&derived ().forest);

      derived ().forest = new_forest;
      derived ().post_adaptation_hook ();
    }
  }

  //=============================================================================
  // Refinement
  //=============================================================================

  /**
   * @brief Harten's neighbour prediction on the leaf level
   *
   * Reference semantics: a significant detail at cell lambda (level L) marks
   * all same-level neighbours of lambda as significant, i.e. their children
   * at L+1 must exist. Leaf formulation used here: every leaf whose family
   * detail is significant (parent in td_set) constructs its same-level face
   * neighbours and pulls the covering leaf up by one level if it is coarser.
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
    auto current_idx = t8_locidx_t {0};
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx) {
      const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
      const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

      for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);
        if (lmi.level () == 0)
          continue;

        // Only leaves of significant families predict their neighbourhood
        if (!derived ().td_set.contains (t8_mra::parent_lmi (lmi)))
          continue;

        const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
        const auto num_faces = scheme->element_get_num_faces (tree_class, element);

        t8_element_t *neigh_element;
        scheme->element_new (tree_class, 1, &neigh_element);

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

        scheme->element_destroy (tree_class, 1, &neigh_element);
      }
    }
  }

  /**
   * @brief One refinement round: analysis + adaptation passes
   *
   * Analysis phase (whole grid, before any adaptation):
   *   1. Compute details of ALL complete leaf families via the
   *      non-destructive two-scale transform (stored at the parent levels).
   *   2. Harten's prediction on each family detail (level L):
   *        - significant (d > eps): same-level neighbours must carry children
   *          -> neighbour_prediction pulls coarser neighbour leaves up
   *        - steep (d > 2^(P+1) eps, L < max_level-1): the family's children
   *          (leaves at L+1) are refined one further level
   * Adaptation phase: one forest pass per level with marked leaves;
   * children data is reconstructed by inverse two-scale with zero details
   * (= exact polynomial subdivision), keeping lmi_map and forest in sync.
   *
   * Leaves at level 0 have no parent family and are never refinement
   * candidates themselves (their families are, via their parents).
   *
   * @param min_level Minimum level to start from
   * @param max_level Maximum level to refine to
   * @return Number of leaves marked for refinement in this round
   */
  unsigned int
  refine_round (int min_level, int max_level)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    derived ().d_map.erase_all ();
    derived ().td_set.erase_all ();
    derived ().refinement_set.erase_all ();
    derived ().coarsening_set.erase_all ();

    //--------------------------------------------------------------------------
    // Analysis phase
    //--------------------------------------------------------------------------

    // 1. Details of all complete leaf families; lmi_map stays untouched.
    derived ().compute_leaf_details (0, max_level);

    // 2. Harten's prediction on the family details
    const auto steep_factor = std::pow (2.0, static_cast<int> (Derived::P_DIM) + 1);
    auto num_families = 0u;
    auto max_ratio = 0.0;

    for (auto L = 0; L < max_level; ++L) {
      for (const auto &[lmi, _] : derived ().d_map[L]) {
        auto detail_norm = derived ().local_detail_norm (lmi);
        for (auto u = 0u; u < Derived::U_DIM; ++u)
          detail_norm[u] /= derived ().c_scaling[u];

        const auto d_max = *std::max_element (detail_norm.begin (), detail_norm.end ());
        const auto local_eps = derived ().c_thresh * derived ().local_threshold_value (lmi);

        ++num_families;
        max_ratio = std::max (max_ratio, d_max / local_eps);

        // Significant: remember for neighbour prediction
        if (d_max > local_eps)
          derived ().td_set.insert (lmi);

        // Steep gradient: refine the family's children (leaves at L+1).
        // Guard keeps the result within max_level.
        if (d_max > steep_factor * local_eps && L < max_level - 1)
          for (const auto &child : t8_mra::children_lmi (lmi))
            derived ().refinement_set.insert (child);
      }
    }

    // 3. Neighbour prediction on the significant families
    neighbour_prediction ();

    // Details are consumed; only the refinement_set is needed from here on
    derived ().d_map.erase_all ();

    auto num_marked = 0u;
    for (auto l = 0; l <= max_level; ++l)
      num_marked += derived ().refinement_set[l].size ();
    t8_debugf ("MRA refine analysis: %u leaf families, %zu significant, max d/eps = %g, %u leaves marked\n",
               num_families, derived ().td_set.size (), max_ratio, num_marked);

    if (num_marked == 0) {
      // Leave no stale state behind (td_set still holds the significant
      // families of the analysis)
      derived ().td_set.erase_all ();
      derived ().refinement_set.erase_all ();
      derived ().coarsening_set.erase_all ();
      return 0;
    }

    //--------------------------------------------------------------------------
    // Adaptation phase: one forest pass per level with marked leaves
    //--------------------------------------------------------------------------

    // Note: starts at min_level (level-0 leaves can be refined; only the
    // analysis needs a parent family, which level 0 does not have)
    for (auto l = min_level; l < max_level; ++l) {
      if (derived ().refinement_set[l].empty ())
        continue;

      t8_debugf ("MRA refine pass %d -> %d: refining %zu elements\n", l, l + 1,
                 derived ().refinement_set[l].size ());

      // Reconstruct children data for the leaves to be refined: inverse
      // two-scale with zero details. This moves the refined leaves' data
      // from level l to their children at level l+1 in lmi_map, exactly
      // mirroring what the forest adapt below does to the leaves.
      for (const auto &lmi : derived ().refinement_set[l])
        derived ().d_map.insert (lmi, element_t {});

      derived ().inverse_multiscale_transformation (l, l + 1);

      // Adapt the forest
      derived ().get_user_data ()->current_refinement_level = l;

      t8_forest_t new_forest;
      t8_forest_ref (derived ().forest);

      new_forest
        = t8_forest_new_adapt (derived ().forest, static_refinement_callback, 0, 0, derived ().get_user_data ());

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
      std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      new_user_data->mra_instance = &derived ();
      new_user_data->current_refinement_level = l;

      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (new_forest, derived ().forest, static_iterate_replace_callback);

      // Cleanup old forest and user data
      auto *old_user_data = derived ().get_user_data ();
      delete old_user_data->lmi_map;
      sc_array_destroy (old_user_data->lmi_idx);
      T8_FREE (old_user_data);
      t8_forest_unref (&derived ().forest);

      derived ().forest = new_forest;
      derived ().post_adaptation_hook ();
    }

    derived ().d_map.erase_all ();
    derived ().td_set.erase_all ();
    derived ().refinement_set.erase_all ();
    derived ().coarsening_set.erase_all ();

    return num_marked;
  }

  /**
   * @brief Perform adaptive refinement from min_level up to max_level
   *
   * Iterates refinement rounds until the grading fixpoint is reached:
   * Harten's neighbour prediction lifts a coarser neighbour of a significant
   * family by one level per round, so a neighbour across a multi-level jump
   * needs several rounds to arrive at the family's leaf level. Newly created
   * children carry zero details and trigger no further marks, hence the
   * iteration terminates (bounded by max_level rounds).
   *
   * @param min_level Minimum level to start from
   * @param max_level Maximum level to refine to
   */
  void
  refine (int min_level, int max_level)
  {
    /// Scaling due to (2.39)
    derived ().c_scaling = derived ().threshold_scaling_factor ();

    for (auto round = 0; round < max_level; ++round) {
      t8_debugf ("MRA refine round %d\n", round);
      if (refine_round (min_level, max_level) == 0)
        break;
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
