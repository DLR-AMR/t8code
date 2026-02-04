/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_forest_adapt.hxx 
 * Definition of a C++ class for the adaptation routine to provide a flexible way of implementing
 * different kinds of adaptation strategies.
 */

#ifndef T8_FOREST_ADAPT_HXX
#define T8_FOREST_ADAPT_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_types.h>

#include <concepts>
#include <functional>
#include <type_traits>

/**
 * Namespace for adaptation related classes and functions.
 */

/* TODO rename to t8_forest_adapt as soon as it is not used as function name anymore.  */
namespace t8_forest_adapt_namespace
{
/** The action to be taken on an element during adaptation. 
   * COARSEN: The element should be coarsened.
   * KEEP: The element should remain as is.
   * REFINE: The element should be refined.
   * Can be extended (e.g., for special refinement types).
  */
class adapt_action {
 public:
  static const int COARSEN = -1;
  static const int KEEP = 0;
  static const int REFINE = 1;
};

/**
   * Callback function type for element adaptation.
   */
using element_callback
  = std::function<adapt_action (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                const t8_scheme *scheme, const t8_eclass_t tree_class)>;

using batched_element_callback
  = std::function<void (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                        const t8_scheme *scheme, const t8_eclass_t tree_class, std::vector<adapt_action> &action)>;

template <typename T>
concept has_element_callback_collect
  = requires (T a, const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions, element_callback cb) {
      {
        a.collect_adapt_actions (forest_from, adapt_actions, cb)
      } -> std::same_as<void>;
    };

template <typename T>
concept has_batched_callback_collect = requires (
  T a, const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions, batched_element_callback cb) {
  {
    a.collect_adapt_actions (forest_from, adapt_actions, cb)
  } -> std::same_as<void>;
};

template <typename T>
concept adapt_actions_collectable = has_element_callback_collect<T> || has_batched_callback_collect<T>;

template <typename T>
concept family_checkable
  = requires (T a, const t8_element_array_t *tree_elements_from, std::vector<const t8_element_t *> &elements_from,
              const t8_locidx_t offset, const t8_scheme *scheme, const t8_eclass_t tree_class) {
      {
        a.family_check (tree_elements_from, elements_from, offset, scheme, tree_class)
      } -> std::convertible_to<bool>;
    };

template <typename T>
concept element_manipulatable = requires (
  T a, t8_element_array_t *elements, const t8_element_array_t *elements_from, const t8_scheme *scheme,
  const t8_eclass_t tree_class, const t8_locidx_t &el_considered, t8_locidx_t &el_inserted, const adapt_action action) {
  {
    a.element_manipulator (elements, elements_from, scheme, tree_class, el_considered, el_inserted, action)
  } -> std::same_as<void>;
};

struct adapt_collector
{
  void
  collect_adapt_actions (const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions,
                         element_callback callback)
  {
    T8_ASSERT (forest_from != nullptr);

    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
    adapt_actions.resize (static_cast<size_t> (local_num_elements));

    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_element_array_t *tree_elements_from = &tree_from->leaf_elements;
      const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (tree_elements_from);
      T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
      const t8_eclass_t tree_class = tree_from->eclass;
      const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

      for (t8_locidx_t el_considered = 0; el_considered < num_el_from; el_considered++) {
        const t8_element_t *element_from = t8_element_array_index_locidx (tree_elements_from, el_considered);
        adapt_actions[el_offset + el_considered] = callback (forest_from, ltree_id, element_from, scheme, tree_class);
      }
      el_offset += num_el_from;
    }
  };

  /** Standard family checker implementation. s*/
  struct family_checker
  {
    bool
    family_check (const t8_element_array_t *tree_elements_from, std::vector<const t8_element_t *> &elements_from,
                  const t8_locidx_t offset, const t8_scheme *scheme, const t8_eclass_t tree_class)
    {
      const int num_siblings = scheme->element_get_num_siblings (tree_class, elements_from[offset]);
      for (int isibling = 0; isibling < num_siblings; isibling++) {
        elements_from[isibling]
          = (const t8_element_t *) t8_element_array_index_locidx (tree_elements_from, offset + (t8_locidx_t) isibling);
        if (scheme->element_get_child_id (tree_class, elements_from[isibling]) != isibling) {
          return false;
        }
      }
      /* elements_are_family expects t8_element_t *const *; build a non-const pointer array */
      std::vector<t8_element_t *> children_nonconst (num_siblings);
      for (int i = 0; i < num_siblings; ++i)
        children_nonconst[i] = const_cast<t8_element_t *> (elements_from[i]);
      const bool is_family = scheme->elements_are_family (tree_class, children_nonconst.data ());
      return is_family;
    }
  };

  struct manipulator
  {
    void
    element_manipulator (t8_element_array_t *elements, const t8_element_array_t *const elements_from,
                         const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t &el_considered,
                         t8_locidx_t &el_inserted, const adapt_action action)
    {
      if (!is_family && action == adapt_action::COARSEN) {
        action = adapt_action::KEEP;
      }
      /* Check that all siblings want to be coarsened */
      if (is_family && action == adapt_action::COARSEN) {
        const auto start = adapt_actions.begin () + static_cast<size_t> (el_offset + el_considered);
        const auto end = start + static_cast<size_t> (num_siblings);
        if (!std::all_of (start, end, [] (const adapt_action &a) { return a == adapt_action::COARSEN; })) {
          action = adapt_action::KEEP;
        }
      }

      switch (action) {
      case adapt_action::COARSEN:
        el_inserted += manipulate_elements<adapt_action::COARSEN> (elements, tree_elements_from, scheme, tree_class,
                                                                   el_inserted, el_offset + el_considered);
        break;
      case adapt_action::KEEP:
        el_inserted += manipulate_elements<adapt_action::KEEP> (elements, tree_elements_from, scheme, tree_class,
                                                                el_inserted, el_offset + el_considered);
        break;
      case adapt_action::REFINE:
        el_inserted += manipulate_elements<adapt_action::REFINE> (elements, tree_elements_from, scheme, tree_class,
                                                                  el_inserted, el_offset + el_considered);
        break;
      default: {
        t8_errorf ("Unknown adapt action.\n");
        SC_ABORT_NOT_REACHED ();
        break;
      }
      }
    }

    /**  * Class implementing a basic adaptation strategy for a forest of trees.
   */
    template <adapt_actions_collectable TCollect, family_checkable TFamily, element_manipulatable TManipulate>
    class adaptor: private TCollect, private TFamily, private TManipulate {
     public:
      /** Constructor for basic_adaptation class.
     * \param [in] forest_in        The forest to be adapted.
     * \param [in] callback_in      The callback function to determine adaptation actions.
     */
      adaptor (t8_forest_t forest, t8_forest_t forest_from, element_callback callback_in)
        : forest (forest), forest_from (forest_from), callback (callback_in)
      {
        T8_ASSERT (forest != nullptr);
        T8_ASSERT (callback);
        if (forest_from != nullptr) {
          t8_forest_ref (forest_from);
        }
        T8_ASSERT (forest != nullptr);
        if (forest != nullptr) {
          t8_forest_ref (forest);
        }
      }

      /** Destructor for adaptor class. */
      ~adaptor ()
      {
        if (forest_from != nullptr) {
          t8_forest_unref (&forest_from);
        }
        if (forest != nullptr) {
          t8_forest_unref (&forest);
        }
      }

      /** Perform the adaptation process on the forest. */
      void
      adapt ();

      /** Type alias for the callback function used in adaptation. 
     * The type depends on whether TCollect uses element_callback or batched_element_callback.
     */
      using callback_type
        = std::conditional_t<has_element_callback_collect<TCollect>, element_callback, batched_element_callback>;

      callback_type callback; /**< The callback function to determine adaptation actions. */
     private:
      /**
     * Profile the adaptation process.
     */
      inline void
      profile_adaptation ()
      {
        T8_ASSERT (forest->profile != nullptr);
        forest->profile->adapt_runtime = -sc_MPI_Wtime ();
      }

      t8_forest_t forest;                      /**< The target forest */
      t8_forest_t forest_from;                 /**< The source forest to adapt from. */
      std::vector<adapt_action> adapt_actions; /**< The adaptation actions for each element in the source forest. */
      bool profiling = false;                  /**< Flag to indicate if profiling is enabled. */
    };
  };

#endif /* T8_FOREST_ADAPT_HXX */
