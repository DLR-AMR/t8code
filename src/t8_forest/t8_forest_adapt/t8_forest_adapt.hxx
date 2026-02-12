/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

namespace t8_adapt
{
/** The action to be taken on an element during adaptation. 
   * COARSEN: The element should be coarsened.
   * KEEP: The element should remain as is.
   * REFINE: The element should be refined.
   * Can be extended (e.g., for special refinement types).
  */
class action {
 public:
  static const int COARSEN = -1; /**< Coarsen the element */
  static const int KEEP = 0;     /**< Keep the element unchanged */
  static const int REFINE = 1;   /**< Refine the element */

  /**
   * Default constructor: KEEP action.
   * \note implicit conversion from int to action.
   */
  action (): value (KEEP)
  {
  }
  /** Constructor from int value.
   * \param [in] v The integer value representing the action.
   * \note implicit conversion from int to action.
   */
  action (int v): value (v)
  {
  }

  /** Conversion operator to int.
   * \return The integer value representing the action.
   * \note implicit conversion from action to int.
   */
  operator int () const
  {
    return value;
  }

  /**
   * Comparison operators with int and action.
   * \param [in] other The other value to compare with.
   * \return True if the values are equal, false otherwise.
   */
  bool
  operator== (int other) const
  {
    return value == other;
  }

  /** Inequality operator with int.
   * \param [in] other The other value to compare with.
   * \return True if the values are not equal, false otherwise.
   */
  bool
  operator!= (int other) const
  {
    return value != other;
  }

  /** Equality operator with another action.
   * \param [in] other The other action to compare with.
   * \return True if the values are equal, false otherwise.
   */
  bool
  operator== (const action &other) const
  {
    return value == other.value;
  }

  /** Inequality operator with another action.
   * \param [in] other The other action to compare with.
   * \return True if the values are not equal, false otherwise.
   */
  bool
  operator!= (const action &other) const
  {
    return value != other.value;
  }

  /**
   * Assignment operator from int.
   * \param [in] v The integer value representing the action to assign.
   * \return A reference to this action after assignment.
   * \note implicit conversion from int to action.  
   */
  action &
  operator= (int v)
  {
    value = v;
    return *this;
  }

 private:
  int value;
};

/** Callback function type for element adaptation. 
 * \param [in] forest        The forest containing the element.
 * \param [in] ltreeid       The local tree ID of the tree containing the element.
 * \param [in] lelement_id   The local element ID of the element to be adapted.
 * \param [in] element       The element to be adapted.
 * \param [in] scheme        The scheme used for the element.
 * \param [in] tree_class    The eclass of the tree containing the element.
 * \return                   The adaptation action to be taken on the element.
*/
using element_callback
  = std::function<action (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_locidx_t lelement_id,
                          const t8_element_t *element, const t8_scheme *scheme, const t8_eclass_t tree_class)>;

/** Callback function type for batched element adaptation. 
 * \param [in] forest        The forest containing the elements.
 * \param [in] ltreeid       The local tree ID of the tree containing the elements.
 * \param [in] element       The array of elements to be adapted.
 * \param [in] scheme        The scheme used for the elements.
 * \param [in] tree_class    The eclass of the tree containing the elements.
 * \param [out] actions      The vector to store the adaptation actions for the elements.
*/
using batched_element_callback
  = std::function<void (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_array_t *element,
                        const t8_scheme *scheme, const t8_eclass_t tree_class, std::vector<action> &actions)>;

/**
 * Concept that detects whether a type TType provides a member function
 * with the signature compatible with:
 *   object.collect_actions(const t8_forest_t, std::vector<action>&, element_callback)
 * returning void.
 *
 * \tparam TType
 *   Type under test. The concept is satisfied when an object `object` of type TType can be
 *   used in an expression
 *     object.collect_actions(forest_from, actions, callback)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - actions is of type std::vector<action>&,
 *     - callback is of type element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename TType>
concept has_element_callback_collect
  = requires (TType object, const t8_forest_t forest_from, std::vector<action> &actions, element_callback callback) {
      {
        object.collect_actions (forest_from, actions, callback)
      } -> std::same_as<void>;
    };

/** Concept that detects whether a type TType provides a member function
 * with the signature compatible with:
 *   object.collect_actions(const t8_forest_t, std::vector<action>&, batched_element_callback)
 * returning void.
 * \tparam TType
 *   Type under test. The concept is satisfied when an object `object` of type TType can be
 *   used in an expression
 *     object.collect_actions(forest_from, actions, callback)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - actions is of type std::vector<action>&,
 *     - callback is of type batched_element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename TType>
concept has_batched_callback_collect = requires (TType object, const t8_forest_t forest_from,
                                                 std::vector<action> &actions, batched_element_callback callback) {
  {
    object.collect_actions (forest_from, actions, callback)
  } -> std::same_as<void>;
};

/** Concept that detects whether a type TType provides a member function
 * with the signature compatible with either:
 *   object.collect_actions(const t8_forest_t, std::vector<action>&, element_callback)
 * or
 *   object.collect_actions(const t8_forest_t, std::vector<action>&, batched_element_callback)
 * returning void.
 * \tparam TType
 *   Type under test. The concept is satisfied when an object `object` of type T can be
 *   used in an expression
 *     object.collect_actions(forest_from, actions, callback)
 *   where:
 *     - forest_from is of type const t8_forest_t,
 *     - actions is of type std::vector<action>&,
 *     - callback is of type element_callback or batched_element_callback,
 *   and the expression is well-formed and yields void.
 */
template <typename TType>
concept actions_collectable = has_element_callback_collect<TType> || has_batched_callback_collect<TType>;

/** Concept that detects whether a type TType provides a member function
 * with the signature compatible with:
 *   object.family_check(const t8_element_array_t*, const t8_locidx_t, const t8_scheme*, const t8_eclass_t)
 * returning bool.
 * \tparam TType
 *   Type under test. The concept is satisfied when an object `object` of type TType can be
 *   used in an expression
 *     object.family_check(tree_elements_from, offset, scheme, tree_class)
 *   where:
 *     - tree_elements_from is of type const t8_element_array_t*,
 *     - offset is of type const t8_locidx_t,
 *     - scheme is of type const t8_scheme*,
 *     - tree_class is of type const t8_eclass_t,
 *   and the expression is well-formed and yields a type convertible to bool.
 */
template <typename TType>
concept family_checkable = requires (TType object, const t8_element_array_t *tree_elements_from,
                                     const t8_locidx_t offset, const t8_scheme *scheme, const t8_eclass_t tree_class) {
  {
    object.family_check (tree_elements_from, offset, scheme, tree_class)
  } -> std::convertible_to<bool>;
};

/** Concept that detects whether a type TType provides a member function
 * with the signature compatible with:
 *  object.element_manipulator(t8_element_array_t*, const t8_element_array_t*, const t8_scheme*,
 *                           const t8_eclass_t, const t8_locidx_t, const t8_locidx_t,
 *                          t8_locidx_t&, const std::vector<action>&, const bool)
 * returning void.
 * \tparam TType
 *   Type under test. The concept is satisfied when an object `object` of type TType can be
 *   used in an expression
 *     object.element_manipulator(elements, elements_from, scheme, tree_class, el_considered, el_offset,
 *                               el_inserted, actions, is_family)
 *   where:
 *     - elements is of type t8_element_array_t*,
 *     - elements_from is of type const t8_element_array_t*,
 *     - scheme is of type const t8_scheme*,
 *     - tree_class is of type const t8_eclass_t,
 *     - el_considered is of type const t8_locidx_t,
 *     - el_offset is of type const t8_locidx_t,
 *     - el_inserted is of type t8_locidx_t&,
 *     - actions is of type const std::vector<action>&,
 *     - is_family is of type const bool,
 *   and the expression is well-formed and yields void.
 */
template <typename TType>
concept element_manipulatable = requires (
  TType object, t8_element_array_t *elements, const t8_element_array_t *const elements_from, const t8_scheme *scheme,
  const t8_eclass_t tree_class, const t8_locidx_t el_considered, const t8_locidx_t el_offset, t8_locidx_t &el_inserted,
  const std::vector<action> &action, const bool is_family) {
  {
    object.element_manipulator (elements, elements_from, scheme, tree_class, el_considered, el_offset, el_inserted,
                                action, is_family)
  } -> std::same_as<void>;
};

/**
   * Policy-based adaptor that drives element-wise adaptation of a target forest from a source forest.
   *
   * This class coordinates the process of collecting adaptation actions from a source forest and applying
   * element-level manipulations to a target forest. It is implemented as a policy-composition class template
   * and privately inherits the provided policy types:
   *  - TCollect: provides collect_actions(forest_from, actions, callback) to produce actions.
   *  - TFamily:  provides family_check(...) to identify family/grouped elements (siblings).
   *  - TManipulate: provides element_manipulator(...) to perform /coarsening/refinement or general element manipulation logic.
   *
   * Template parameters
   * \tparam TCollect    A type satisfying actions_collectable: must expose collect_actions.
   * \tparam TFamily     A type satisfying family_checkable: must expose family_check to detect families.
   * \tparam TManipulate A type satisfying element_manipulatable: must expose element_manipulator.
   * \tparam TRecursive  A boolean flag (default false) that can be used by policies to enable recursive adaptation logic.
   *
   * Overview
   * The adaptor holds references to a "target" forest and a "source" forest (forest and forest_from). On construction it
   * retains references (increments reference counts) for any non-null forest handles; on destruction it releases them.
   * The adapt() method:
   *  - Asserts valid state and optionally starts profiling.
   *  - Uses the TCollect policy to populate actions for each element in the source forest.
   *  - Iterates over local trees of the source forest, and for each tree:
   *      - Retrieves corresponding tree data from both source and target forests.
   *      - Accesses the source tree's leaf element array and the number of source elements.
   *      - For each considered source element (taking element siblings/families into account), uses the TFamily
   *        policy to determine whether the current elements form a family, reads the precomputed adapt action for
   *        the element, and calls the TManipulate policy to update the target tree's element array accordingly.
   *      - Maintains per-tree and global offsets so that actions are applied using a linear index across
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
   * - inline void profile_adaptation_start()
   *     - Starts adaptation timing by writing into forest->profile->adapt_runtime using the MPI wall-time helper.
   *     - Assumes a non-null profiling structure on the forest.
   * - inline void profile_adaptation_end()
   *    - Ends adaptation timing by computing the elapsed time since profile_adaptation_start and accumulating it
   *      into forest->profile->adapt_runtime.
   * - void recursive_adapt(...)
   *    - A helper function that can be used when TRecursive is true to perform recursive adaptation logic. It takes parameters describing the current element being considered,
   *      the source and target element arrays, the scheme and tree class, the current offset and count of inserted elements, the precomputed actions, and whether the current element is part of a family. 
   *      It manages the recursive insertion of child elements via a depth-first traversal, applying the TManipulate policy at each level and using the callback to determine actions for child elements. 
   *      This function is only enabled when TRecursive is true and TCollect uses a non-batched element_callback.
   *
   * Data members
   * - t8_forest_t forest
   *     - The target forest that will be modified by adaptation operations.
   * - t8_forest_t forest_from
   *     - The source forest from which adaptation decisions are derived.
   * - std::vector<action> actions
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
template <actions_collectable TCollect, family_checkable TFamily, element_manipulatable TManipulate,
          bool TRecursive = false>
class adaptor: private TCollect, private TFamily, private TManipulate {
 public:
  /** The type of callback used for collecting adaptation actions. */
  using callback_type
    = std::conditional_t<has_element_callback_collect<TCollect>, element_callback, batched_element_callback>;

  /** Constructor for basic_adaptation class.
   * \param [in] forest            The target forest to be adapted.
   * \param [in] forest_from       The source forest to adapt from.
   * \param [in] callback_in       The callback function to determine adaptation actions for elements.
   * \param [in] profiling_in      Flag to indicate if profiling should be enabled during adaptation.
   *
   * \note The constructor increments reference counts for non-null forest handles, and the destructor will release them.
   */
  adaptor (t8_forest_t forest, t8_forest_t forest_from, callback_type callback_in, bool profiling_in = false)
    : callback (callback_in), forest (forest), forest_from (forest_from), profiling (profiling_in)
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
      profile_adaptation_start ();
    }
    T8_ASSERT (forest_from != nullptr);

    TCollect::collect_actions (forest_from, actions, callback);

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
      if (num_el_from > 0) {
        /* index of the elements in source tree */
        t8_locidx_t el_considered = 0;
        while (el_considered < num_el_from) {
          t8_locidx_t el_inserted = 0;
          const bool is_family = TFamily::family_check (tree_elements_from, el_considered, scheme, tree_class);
          if constexpr (TRecursive) {
            t8_debugf ("[D] recursive adapt: tree %d, element %d, is_family %d\n", ltree_id, el_considered, is_family);
            recursive_adapt (elements, tree_elements_from, scheme, tree_class, ltree_id, el_considered, el_offset,
                             el_inserted, actions, is_family);
          }
          else {
            t8_debugf ("[D] non-recursive adapt: tree %d, element %d, is_family %d\n", ltree_id, el_considered,
                       is_family);
            /* manipulator step*/
            TManipulate::element_manipulator (elements, tree_elements_from, scheme, tree_class, el_considered,
                                              el_offset, el_inserted, actions, is_family);
          }
          el_considered++;
          el_offset += el_inserted;
          forest->local_num_leaf_elements += el_inserted;
        }
      }
      tree->elements_offset = el_offset;
    }
    if (profiling) {
      profile_adaptation_end ();
    }
  }

  callback_type callback; /**< The callback function to determine adaptation actions. */
 private:
  /* WARNING: Currently only available with iterative callback*/
  void
  recursive_adapt (t8_element_array_t *elements, const t8_element_array_t *elements_from, const t8_scheme *scheme,
                   const t8_eclass_t tree_class, const t8_locidx_t ltree_id, const t8_locidx_t el_considered,
                   const t8_locidx_t el_offset, t8_locidx_t &el_inserted, const std::vector<action> &actions,
                   const bool is_family)
    requires (TRecursive && has_element_callback_collect<TCollect>)
  {
    const size_t element_size = scheme->get_element_size (tree_class);
    t8_locidx_t recursively_inserted = 0;
    TManipulate::element_manipulator (elements, elements_from, scheme, tree_class, el_considered, el_offset,
                                      recursively_inserted, actions, is_family);
    t8_element_array_t *children = t8_element_array_new (scheme, tree_class);
    /* add all elements from el_offset to el_offset + recursively_inserted */
    for (t8_locidx_t i = 0; i < recursively_inserted; ++i) {
      t8_element_t *child = t8_element_array_push (children);
      child = t8_element_array_index_locidx_mutable (elements, el_offset + i);
      memcpy ((void *) child, (void *) t8_element_array_index_locidx (elements_from, el_offset + el_considered),
              element_size);
    }
    std::vector<action> child_actions;
    while (t8_element_array_get_count (children) > 0) {
      const int next_child = t8_element_array_get_count (children) - 1;
      const t8_element_t *child = t8_element_array_pop (children);
      t8_debugf ("[D] recursive adapt: tree %d, child element %d\n", ltree_id, next_child);
      /** TODO: next_child is currently the wrong linear id of child */
      child_actions.push_back (callback (forest_from, ltree_id, next_child, child, scheme, tree_class));
      const int level = scheme->element_get_level (tree_class, child);
      if (child_actions[next_child] != action::COARSEN && level < forest->maxlevel) {
        t8_locidx_t inserted = recursively_inserted;
        TManipulate::element_manipulator (elements, children, scheme, tree_class, next_child, el_offset,
                                          recursively_inserted, child_actions, false);
        inserted = recursively_inserted - inserted;
        /* push children */
        t8_element_array_push_count (children, inserted);
        for (t8_locidx_t i = recursively_inserted - 1; i >= recursively_inserted - inserted; --i) {
          t8_element_t *new_child = t8_element_array_index_locidx_mutable (children, i);
          memcpy ((void *) new_child, (void *) t8_element_array_index_locidx (elements, el_offset + i), element_size);
        }
      }
      if (child_actions[next_child] != action::REFINE) {
        child_actions.pop_back ();
      }
    }

    el_inserted += recursively_inserted;
  }

  /**
   * Profile the adaptation process.
   */
  inline void
  profile_adaptation_start ()
  {
    T8_ASSERT (forest->profile != nullptr);
    forest->profile->adapt_runtime = -sc_MPI_Wtime ();
  }

  /**
   * End profiling the adaptation process.
   */
  inline void
  profile_adaptation_end ()
  {
    T8_ASSERT (forest->profile != nullptr);
    forest->profile->adapt_runtime += sc_MPI_Wtime ();
  }

  t8_forest_t forest;          /**< The target forest */
  t8_forest_t forest_from;     /**< The source forest to adapt from. */
  std::vector<action> actions; /**< The adaptation actions for each element in the source forest. */
  bool profiling = false;      /**< Flag to indicate if profiling is enabled. */
};                             // class adaptor

};     // namespace t8_adapt
#endif /* T8_FOREST_ADAPT_HXX */
