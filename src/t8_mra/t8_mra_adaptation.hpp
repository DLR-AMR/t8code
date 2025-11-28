#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/t8_mra_base.hpp"
#include "t8_forest/t8_forest_adapt.h"
#include "t8_forest/t8_forest_ghost.h"
#include "t8_forest/t8_forest_iterate.h"

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
    const auto offset = t8_forest_get_tree_element_offset (forest, which_tree);
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
    using levelmultiindex = typename Derived::levelmultiindex;

    auto *user_data = reinterpret_cast<t8_mra::forest_data<element_t> *> (t8_forest_get_user_data (forest_from));

    // Only process elements at current refinement level
    const auto element_level = scheme->element_get_level (tree_class, elements[0]);
    if (element_level != user_data->current_refinement_level)
      return 0;

    // Check if element is marked for refinement
    const auto tree_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest_from, which_tree);
    const auto lmi = levelmultiindex (tree_id, elements[0], scheme);

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
      // Refinement: set children LMIs
      const auto children = t8_mra::children_lmi (old_lmi);
      for (int i = 0; i < num_incoming; ++i)
        t8_mra::set_lmi_forest_data (new_user_data, first_incoming + i, children[i]);
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
   * ONE-TO-ONE implementation of old t8_mra.hpp::coarsening_new()
   *
   * @param min_level Minimum level to coarsen to
   * @param max_level Maximum level to start from
   */
  void
  coarsening_new (int min_level, int max_level)
  {
    using element_t = typename Derived::element_t;
    using levelmultiindex = typename Derived::levelmultiindex;

    // Static lambda callbacks (must be static to capture 'this')
    static auto static_coarsening_callback
      = [this] (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                t8_locidx_t local_ele_idx, const t8_scheme_c *scheme, const int is_family, const int num_elements,
                t8_element_t *elements[]) -> int {
      return static_cast<Derived *> (this)->coarsening_callback (
        forest, forest_from, which_tree, tree_class, local_ele_idx, scheme, is_family, num_elements, elements);
    };

    static auto static_iterate_replace_callback
      = [this] (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                const t8_scheme *scheme, int refine, int num_outgoing, t8_locidx_t first_outgoing, int num_incoming,
                t8_locidx_t first_incoming) -> void {
      return static_cast<Derived *> (this)->iterate_replace_callback (forest_old, forest_new, which_tree, tree_class,
                                                                      scheme, refine, num_outgoing, first_outgoing,
                                                                      num_incoming, first_incoming);
    };

    /// Scaling due to (2.39)
    derived ().c_scaling = derived ().threshold_scaling_factor ();

    for (auto l = max_level; l > min_level; --l) {
      t8_forest_t new_forest;
      t8_forest_ref (derived ().forest);

      derived ().get_user_data ()->current_refinement_level = l;

      std::cout << "Before mst: " << derived ().get_user_data ()->lmi_map->size () << "\n";
      derived ().multiscale_transformation (l - 1, l);
      derived ().threshold (l - 1, l);
      // restore_balancing (l - 1, l); /// TODO
      derived ().generate_td_tree (l - 1, l);
      derived ().sync_d_with_td (min_level, max_level);

      derived ().inverse_multiscale_transformation (l - 1, l);

      std::cout << "After mst: " << derived ().get_user_data ()->lmi_map->size () << "\n";
      std::cout << "  Level " << (l - 1)
                << " in map: " << derived ().get_user_data ()->lmi_map->operator[] (l - 1).size () << "\n";
      std::cout << "  Level " << l << " in map: " << derived ().get_user_data ()->lmi_map->operator[] (l).size ()
                << "\n";
      std::cout << "  coarsening_set[" << l << "] size: " << derived ().coarsening_set[l].size () << "\n";

      new_forest = t8_forest_new_adapt (
        derived ().forest,
        [] (auto *forest, auto *forest_from, auto which_tree, auto tree_class, auto local_ele_idx, auto *scheme,
            const auto is_family, const auto num_elements, auto *elements[]) -> int {
          return static_coarsening_callback (forest, forest_from, which_tree, tree_class, local_ele_idx, scheme,
                                             is_family, num_elements, elements);
        },
        0, 0, derived ().get_user_data ());

      t8_mra::forest_data<element_t> *new_user_data;
      new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);

      new_user_data->lmi_map = new t8_mra::levelindex_map<levelmultiindex, element_t> (derived ().maximum_level);
      std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

      const auto num_new_local_elements = t8_forest_get_local_num_leaf_elements (new_forest);
      const auto num_new_ghost_elements = t8_forest_get_num_ghosts (new_forest);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (levelmultiindex), num_new_local_elements + num_new_ghost_elements);
      t8_forest_set_user_data (new_forest, new_user_data);
      t8_forest_iterate_replace (
        new_forest, derived ().forest,
        [] (auto *forest_old, auto *forest_new, auto which_tree, const auto tree_class, const auto *scheme, auto refine,
            auto num_outgoing, auto first_outgoing, auto num_incoming, t8_locidx_t first_incoming) -> void {
          static_iterate_replace_callback (forest_old, forest_new, which_tree, tree_class, scheme, refine, num_outgoing,
                                           first_outgoing, num_incoming, first_incoming);
        });

      // Debug: Count elements at each level in the forest
      // {
      //   std::map<int, int> level_counts;
      //   auto *new_forest_data = t8_mra::get_mra_forest_data<element_t> (new_forest);
      //   const auto *new_scheme = t8_forest_get_scheme (new_forest);
      //
      //   for (t8_locidx_t tree_idx = 0; tree_idx < t8_forest_get_num_local_trees (new_forest); ++tree_idx) {
      //     const auto tree_class = t8_forest_get_tree_class (new_forest, tree_idx);
      //     const auto num_elems = t8_forest_get_tree_num_leaf_elements (new_forest, tree_idx);
      //     for (t8_locidx_t elem_idx = 0; elem_idx < num_elems; ++elem_idx) {
      //       const auto *elem = t8_forest_get_leaf_element_in_tree (new_forest, tree_idx, elem_idx);
      //       const int elem_level = new_scheme->element_get_level (tree_class, elem);
      //       level_counts[elem_level]++;
      //     }
      //   }
      //
      //   std::cout << "  Forest elements by level: ";
      //   for (const auto &[lev, cnt] : level_counts) {
      //     std::cout << "L" << lev << "=" << cnt << " ";
      //   }
      //   std::cout << "\n";
      // }

      derived ().d_map.erase_all ();
      derived ().td_set.erase_all ();
      derived ().refinement_set.erase_all ();
      derived ().coarsening_set.erase_all ();

      // cleanup() is called BEFORE updating forest pointer (old version does this)
      // Note: In the new architecture, we manually clean up old_user_data instead
      auto *old_user_data = derived ().get_user_data ();
      delete old_user_data->lmi_map;
      sc_array_destroy (old_user_data->lmi_idx);
      T8_FREE (old_user_data);
      t8_forest_unref (&derived ().forest);

      derived ().forest = new_forest;

      // Update vertex orders after adaptation
      derived ().post_adaptation_hook ();
    }
  }

  // void
  // coarsening_new (int min_level, int max_level)
  // {
  //   using element_t = typename Derived::element_t;
  //
  //   // Compute scaling factors (2.39)
  //   derived ().c_scaling = derived ().threshold_scaling_factor ();
  //
  //   for (auto l = max_level; l > min_level; --l) {
  //     // Set current level for callbacks
  //     derived ().get_user_data ()->current_refinement_level = l;
  //
  //     // Forward MST: compute parent coefficients and details
  //     derived ().multiscale_transformation (l - 1, l);
  //
  //     // Threshold details to determine which elements to coarsen
  //     derived ().threshold (l - 1, l);
  //
  //     // Generate td_tree: ensure parents of significant elements are included
  //     derived ().generate_td_tree (l - 1, l);
  //
  //     // Sync d_map with td_set: remove insignificant parents, mark children for coarsening
  //     derived ().sync_d_with_td (min_level, max_level);
  //
  //     // Inverse MST: reconstruct children for significant parents (they stay refined)
  //     derived ().inverse_multiscale_transformation (l - 1, l);
  //
  //     std::cout << "Level " << l << " -> " << (l - 1) << ": "
  //               << "Details: " << derived ().d_map[l - 1].size () << " parents, "
  //               << "Coarsening: " << derived ().coarsening_set[l].size () << " children\n";
  //
  //     // Create adapted forest using t8_forest_new_adapt
  //     t8_forest_ref (derived ().forest);
  //     t8_forest_t forest_new = t8_forest_new_adapt (derived ().forest, static_coarsening_callback, 0, 0,
  //                                                     derived ().get_user_data ());
  //
  //     // Apply balancing if requested
  //     if (derived ().balanced) {
  //       t8_forest_t balanced_forest;
  //       t8_forest_init (&balanced_forest);
  //       t8_forest_set_balance (balanced_forest, forest_new, 0);
  //       t8_forest_commit (balanced_forest);
  //       forest_new = balanced_forest;
  //     }
  //
  //     // Allocate new user data and swap lmi_map
  //     auto *new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
  //     new_user_data->lmi_map = new t8_mra::levelindex_map<typename Derived::levelmultiindex, element_t> (
  //       derived ().maximum_level);
  //     std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);
  //
  //     // Allocate lmi_idx array
  //     const auto num_elements = t8_forest_get_local_num_leaf_elements (forest_new);
  //     const auto num_ghosts = t8_forest_get_num_ghosts (forest_new);
  //     new_user_data->lmi_idx = sc_array_new_count (sizeof (typename Derived::levelmultiindex),
  //                                                   num_elements + num_ghosts);
  //     new_user_data->mra_instance = &derived ();
  //     new_user_data->current_refinement_level = l;
  //
  //     t8_forest_set_user_data (forest_new, new_user_data);
  //
  //     // Copy element data from old to new forest
  //     t8_forest_iterate_replace (forest_new, derived ().forest, static_iterate_replace_callback);
  //
  //     // Cleanup old forest and user data
  //     auto *old_user_data = derived ().get_user_data ();
  //     delete old_user_data->lmi_map;
  //     sc_array_destroy (old_user_data->lmi_idx);
  //     T8_FREE (old_user_data);
  //     t8_forest_unref (&derived ().forest);
  //
  //     // Update forest pointer
  //     derived ().forest = forest_new;
  //
  //     // Call post-adaptation hook (for element-specific cleanup like vertex order updates)
  //     derived ().post_adaptation_hook ();
  //
  //     // Clear sets for next iteration
  //     derived ().d_map.erase_all ();
  //     derived ().td_set.erase_all ();
  //     derived ().refinement_set.erase_all ();
  //     derived ().coarsening_set.erase_all ();
  //   }
  // }

  //=============================================================================
  // Refinement
  //=============================================================================

  /**
   * @brief Perform adaptive refinement from min_level up to max_level
   *
   * @param min_level Minimum level to start from
   * @param max_level Maximum level to refine to
   */
  void
  refinement_new (int min_level, int max_level)
  {
    using element_t = typename Derived::element_t;

    for (auto l = min_level; l < max_level; ++l) {
      // Set current level for callbacks
      derived ().get_user_data ()->current_refinement_level = l;

      // Sync detail map with significant details set
      derived ().sync_d_with_td (l, l + 1);

      std::cout << "Level " << l << " -> " << (l + 1) << ": "
                << "Refining: " << derived ().refinement_set[l].size () << " elements\n";

      // Inverse MST: reconstruct children for elements to be refined
      derived ().inverse_multiscale_transformation (l, l + 1);

      // Create adapted forest using t8_forest_new_adapt
      t8_forest_ref (derived ().forest);
      t8_forest_t forest_new
        = t8_forest_new_adapt (derived ().forest, static_refinement_callback, 0, 0, derived ().get_user_data ());

      // Apply balancing if requested
      if (derived ().balanced) {
        t8_forest_t balanced_forest;
        t8_forest_init (&balanced_forest);
        t8_forest_set_balance (balanced_forest, forest_new, 0);
        t8_forest_commit (balanced_forest);
        forest_new = balanced_forest;
      }

      // Allocate new user data and swap lmi_map
      auto *new_user_data = T8_ALLOC (t8_mra::forest_data<element_t>, 1);
      new_user_data->lmi_map
        = new t8_mra::levelindex_map<typename Derived::levelmultiindex, element_t> (derived ().maximum_level);
      std::swap (new_user_data->lmi_map, derived ().get_user_data ()->lmi_map);

      // Allocate lmi_idx array
      const auto num_elements = t8_forest_get_local_num_leaf_elements (forest_new);
      const auto num_ghosts = t8_forest_get_num_ghosts (forest_new);
      new_user_data->lmi_idx
        = sc_array_new_count (sizeof (typename Derived::levelmultiindex), num_elements + num_ghosts);
      new_user_data->mra_instance = &derived ();
      new_user_data->current_refinement_level = l;

      t8_forest_set_user_data (forest_new, new_user_data);

      // Copy element data from old to new forest
      t8_forest_iterate_replace (forest_new, derived ().forest, static_iterate_replace_callback);

      // Cleanup old forest and user data
      auto *old_user_data = derived ().get_user_data ();
      delete old_user_data->lmi_map;
      sc_array_destroy (old_user_data->lmi_idx);
      T8_FREE (old_user_data);
      t8_forest_unref (&derived ().forest);

      // Update forest pointer
      derived ().forest = forest_new;

      // Call post-adaptation hook (for element-specific cleanup like vertex order updates)
      derived ().post_adaptation_hook ();

      // Clear sets for next iteration
      derived ().d_map.erase_all ();
      derived ().td_set.erase_all ();
      derived ().refinement_set.erase_all ();
      derived ().coarsening_set.erase_all ();
    }
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
