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

/** Callback function type for element adaptation. 
 * \param [in] forest        The forest containing the element.
 * \param [in] ltreeid       The local tree ID of the tree containing the element.
 * \param [in] element       The element to be adapted.
 * \param [in] scheme        The scheme used for the element.
 * \param [in] tree_class    The eclass of the tree containing the element.
 * \return                   The adaptation action to be taken on the element.
*/
using element_callback
  = std::function<adapt_action (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                const t8_scheme *scheme, const t8_eclass_t tree_class)>;

/** Callback function type for batched element adaptation. 
 * \param [in] forest        The forest containing the elements.
 * \param [in] ltreeid       The local tree ID of the tree containing the elements.
 * \param [in] element       The array of elements to be adapted.
 * \param [in] scheme        The scheme used for the elements.
 * \param [in] tree_class    The eclass of the tree containing the elements.
 * \param [out] action       The vector to store the adaptation actions for the elements.
*/
using batched_element_callback
  = std::function<void (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_array_t *element,
                        const t8_scheme *scheme, const t8_eclass_t tree_class, std::vector<adapt_action> &action)>;

/**
 * Concept that detects whether a type T provides a member function
 * with the signature compatible with:
 *   a.collect_adapt_actions(const t8_forest_t, std::vector<adapt_action>&, element_callback)
 * returning void.
 *
 * @tparam T
 *   Type under test. The concept is satisfied when an object `a` of type T can be
 *   used in an expression
 *     a.collect_adapt_actions(forest_from, adapt_actions, cb)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - adapt_actions is of type std::vector<adapt_action>&,
 *     - cb is of type element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename T>
concept has_element_callback_collect
  = requires (T a, const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions, element_callback cb) {
      {
        a.collect_adapt_actions (forest_from, adapt_actions, cb)
      } -> std::same_as<void>;
    };

/** Concept that detects whether a type T provides a member function
 * with the signature compatible with:
 *   a.collect_adapt_actions(const t8_forest_t, std::vector<adapt_action>&, batched_element_callback)
 * returning void.
 * @tparam T
 *   Type under test. The concept is satisfied when an object `a` of type T can be
 *   used in an expression
 *     a.collect_adapt_actions(forest_from, adapt_actions, cb)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - adapt_actions is of type std::vector<adapt_action>&,
 *     - cb is of type batched_element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename T>
concept has_batched_callback_collect = requires (
  T a, const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions, batched_element_callback cb) {
  {
    a.collect_adapt_actions (forest_from, adapt_actions, cb)
  } -> std::same_as<void>;
};

/** Concept that detects whether a type T provides a member function
 * with the signature compatible with either:
 *   a.collect_adapt_actions(const t8_forest_t, std::vector<adapt_action>&, element_callback)
 * or
 *   a.collect_adapt_actions(const t8_forest_t, std::vector<adapt_action>&, batched_element_callback)
 * returning void.
 * @tparam T
 *   Type under test. The concept is satisfied when an object `a` of type T can be
 *   used in an expression
 *     a.collect_adapt_actions(forest_from, adapt_actions, cb)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - adapt_actions is of type std::vector<adapt_action>&,
 *     - cb is of type element_callback or batched_element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename T>
concept adapt_actions_collectable = has_element_callback_collect<T> || has_batched_callback_collect<T>;

/** Concept that detects whether a type T provides a member function
 * with the signature compatible with:
 *   a.family_check(const t8_element_array_t*, std::vector<const t8_element_t*>&, const t8_locidx_t,
 *                  const t8_scheme*, const t8_eclass_t)
 * returning bool.
 * @tparam T
 *   Type under test. The concept is satisfied when an object `a` of type T can be
 *   used in an expression
 *     a.family_check(tree_elements_from, elements_from, offset, scheme, tree_class)
 *   where:
 *     - tree_elements_from is of type const t8_element_array_t*,
 *     - elements_from is of type std::vector<const t8_element_t*>&,
 *     - offset is of type const t8_locidx_t,
 *     - scheme is of type const t8_scheme*,
 *     - tree_class is of type const t8_eclass_t,
 *   and the expression is well-formed and yields a type convertible to bool.
 */
template <typename T>
concept family_checkable
  = requires (T a, const t8_element_array_t *tree_elements_from, std::vector<const t8_element_t *> &elements_from,
              const t8_locidx_t offset, const t8_scheme *scheme, const t8_eclass_t tree_class) {
      {
        a.family_check (tree_elements_from, elements_from, offset, scheme, tree_class)
      } -> std::convertible_to<bool>;
    };

/** Concept that detects whether a type T provides a member function
 * with the signature compatible with:
 *   a.element_manipulator(t8_element_array_t*, const t8_element_array_t*,
 *                         const t8_scheme*, const t8_eclass_t, const t8_locidx_t&,
 *                         t8_locidx_t&, const adapt_action)
 * returning void.
 * @tparam T
 *   Type under test. The concept is satisfied when an object `a` of type T can be
 *   used in an expression
 *     a.element_manipulator(elements, elements_from, scheme, tree_class, el_considered, el_inserted, action)
 *   where:
 *     - elements is of type t8_element_array_t*,
 *     - elements_from is of type const t8_element_array_t*,
 *     - scheme is of type const t8_scheme*,
 *     - tree_class is of type const t8_eclass_t,
 *     - el_considered is of type const t8_locidx_t&,
 *     - el_inserted is of type t8_locidx_t&,
 *     - action is of type const adapt_action,
 *   and the expression is well-formed and yields void.
 */
template <typename T>
concept element_manipulatable = requires (
  T a, t8_element_array_t *elements, const t8_element_array_t *elements_from, const t8_scheme *scheme,
  const t8_eclass_t tree_class, const t8_locidx_t &el_considered, t8_locidx_t &el_inserted, const adapt_action action) {
  {
    a.element_manipulator (elements, elements_from, scheme, tree_class, el_considered, el_inserted, action)
  } -> std::same_as<void>;
};

/** Standard adapt action collector implementation. */
struct adapt_collector
{
  /** Collect adapt actions for all elements in the forest.
   * \param [in] forest_from       The forest containing the elements to be adapted.
   * \param [out] adapt_actions    The vector to store the adapt actions for all elements.
   * \param [in] callback          The callback function to determine the adapt action for each element.
   */
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
  }
};

/** Standard family checker implementation. */
struct family_checker
{
  /**
     * Check if the elements in the given array are siblings0 in the tree and form a family.
     * \param [in] tree_elements_from The array of elements in the tree.
     * \param [out] elements_from      The vector to store the elements from the tree.
     * \param [in] offset              The offset to start checking from.
     * \param [in] scheme              The scheme to use for checking.
     * \param [in] tree_class          The class of the tree.
     * \return True if the elements are siblings and form a family, false otherwise.
     */
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

/**
   * Class implementing a basic element manipulation strategy.
   */
struct manipulator
{
  /** Manipulate elements based on the given adapt action.
     * \param [in,out] elements          The array of elements to be manipulated.
     * \param [in] tree_elements_from    The array of elements from the source tree.
     * \param [in] scheme                The scheme to use for manipulation.
     * \param [in] tree_class            The class of the tree.
     * \param [in] el_considered         The index of the element being considered.
     * \param [in,out] el_inserted       The index of the next element to be inserted.
     * \param [in] action                The adapt action to be performed.
     */
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
  };
};

/**
   * Policy-based adaptor that drives element-wise adaptation of a target forest from a source forest.
   *
   * This class coordinates the process of collecting adaptation actions from a source forest and applying
   * element-level manipulations to a target forest. It is implemented as a policy-composition class template
   * and privately inherits the provided policy types:
   *  - TCollect: provides collect_adapt_actions(forest_from, adapt_actions, callback) to produce actions.
   *  - TFamily:  provides family_check(...) to identify family/grouped elements (siblings).
   *  - TManipulate: provides element_manipulator(...) to perform /coarsening/refinement or general element manipulation logic.
   *
   * Template parameters
   * @tparam TCollect    A type satisfying adapt_actions_collectable: must expose collect_adapt_actions.
   * @tparam TFamily     A type satisfying family_checkable: must expose family_check to detect families.
   * @tparam TManipulate A type satisfying element_manipulatable: must expose element_manipulator.
   *
   * Overview
   * The adaptor holds references to a "target" forest and a "source" forest (forest and forest_from). On construction it
   * retains references (increments reference counts) for any non-null forest handles; on destruction it releases them.
   * The adapt() method:
   *  - Asserts valid state and optionally starts profiling.
   *  - Uses the TCollect policy to populate adapt_actions for each element in the source forest.
   *  - Iterates over local trees of the source forest, and for each tree:
   *      - Retrieves corresponding tree data from both source and target forests.
   *      - Accesses the source tree's leaf element array and the number of source elements.
   *      - For each considered source element (taking element siblings/families into account), uses the TFamily
   *        policy to determine whether the current elements form a family, reads the precomputed adapt action for
   *        the element, and calls the TManipulate policy to update the target tree's element array accordingly.
   *      - Maintains per-tree and global offsets so that adapt_actions are applied using a linear index across
   *        all source elements.
   *
   * Ownership and lifetime
   * - The adaptor stores raw t8_forest_t handles for both forest and forest_from. When constructed it will call
   *   t8_forest_ref on non-null inputs and t8_forest_unref on destruction. Callers should treat the adaptor as
   *   owning an additional reference to the provided forest handles for the lifetime of the adaptor.
   *
   * Public interface
   * - adaptor(forest, forest_from, callback_in)
   *     - Constructs an adaptor for adapting `forest` from `forest_from` using `callback_in` to decide actions.
   *     - Preconditions: `callback_in` must be valid; `forest` is asserted non-null.
   *     - Increments reference counts for non-null forest handles.
   *
   * - ~adaptor()
   *     - Releases retained references to the forests.
   *
   * - void adapt()
   *     - Performs the full adaptation: collects actions via TCollect, then visits source elements and applies
   *       element-level manipulations via TManipulate. Updates per-tree element offsets in the target forest.
   *     - If profiling is enabled, profile_adaptation() is invoked to start timing.
   *
   * - using callback_type = std::conditional_t<has_element_callback_collect<TCollect>, element_callback, batched_element_callback>
   *     - Alias describing the type of callback expected by the collect policy. The stored member `callback` uses
   *       this alias so the adaptor can support either single-element or batched callbacks depending on TCollect.
   *
   * Protected / private helpers
   * - inline void profile_adaptation()
   *     - Starts adaptation timing by writing into forest->profile->adapt_runtime using the MPI wall-time helper.
   *     - Assumes a non-null profiling structure on the forest.
   *
   * Data members
   * - t8_forest_t forest
   *     - The target forest that will be modified by adaptation operations.
   * - t8_forest_t forest_from
   *     - The source forest from which adaptation decisions are derived.
   * - std::vector<adapt_action> adapt_actions
   *     - Linear array of adaptation actions for each element in the source forest. Populated by the TCollect policy.
   * - callback_type callback
   *     - The callback provided by the caller; forwarded to the collect policy.
   * - bool profiling
   *     - When true, the adaptor calls profile_adaptation() before starting adaptation.
   *
   * Notes and implementation considerations
   * - The adaptor is designed to separate responsibilities: collection of adapt decisions, detection of sibling/family
   *   groupings, and the low-level manipulation of element arrays are all delegated to policy types. This allows
   *   different collection and manipulation strategies to be plugged in without changing the control flow.
   * - The adaptor relies on the forest and tree data structures exposing stable array indexing via t8_element_array_* APIs.
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
  adapt ()
  {
    T8_ASSERT (forest != nullptr);
    if (profiling) {
      profile_adaptation ();
    }
    T8_ASSERT (forest_from != nullptr);

    TCollect::collect_adapt_actions (forest_from, adapt_actions, callback);

    /* Offset per tree in the source forest */
    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    /* Get the scheme used by the forest */
    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      /* get the trees from both forests. */
      t8_tree_t tree = t8_forest_get_tree (forest, ltree_id);
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      /* get the leaf arrays from both forests */
      t8_element_array_t *elements = &tree->leaf_elements;
      const t8_element_array_t *tree_elements_from = &tree_from->leaf_elements;
      /* Get the number of elements in the source tree */
      const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (tree_elements_from);
      T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
      const t8_eclass_t tree_class = tree_from->eclass;
      /* Continue only if tree_from is not empty */
      if (num_el_from < 0) {
        const t8_element_t *first_element_from = t8_element_array_index_locidx (tree_elements_from, 0);
        t8_locidx_t curr_size_elements_from = scheme->element_get_num_siblings (tree_class, first_element_from);
        /* index of the elements in source tree */
        t8_locidx_t el_considered = 0;
        /* index of the elements in target tree */
        t8_locidx_t el_inserted = 0;
        std::vector<const t8_element_t *> elements_temp;

        while (el_considered < num_el_from) {
          const t8_locidx_t num_siblings = scheme->element_get_num_siblings (
            tree_class, t8_element_array_index_locidx (tree_elements_from, el_considered));
          if (num_siblings > curr_size_elements_from) {
            elements_temp.resize (num_siblings);
            curr_size_elements_from = num_siblings;
          }
          for (int isibling = 0; isibling < num_siblings && el_considered + isibling < num_el_from; isibling++) {
            elements_temp[isibling] = (const t8_element_t *) t8_element_array_index_locidx (
              tree_elements_from, el_considered + (t8_locidx_t) isibling);
            if (scheme->element_get_child_id (tree_class, elements_temp[isibling]) != isibling) {
              break;
            }
          }

          const bool is_family
            = TFamily::family_check (tree_elements_from, elements_temp, el_considered, scheme, tree_class);
          const adapt_action action = adapt_actions[el_offset + el_considered];

          /* manipulator step*/
          TManipulate::element_manipulator (elements, tree_elements_from, scheme, tree_class, el_considered,
                                            el_inserted, action, is_family);
          el_considered++;
        }
      }
      tree->elements_offset = el_offset;
      el_offset += num_el_from;
    }
  }
  /** The type of callback used for collecting adaptation actions. */
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
}
;
#endif /* T8_FOREST_ADAPT_HXX */
