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

#pragma once

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>
#include <vector>

namespace t8_standard_adapt
{
/** Standard adapt action collector implementation. */
struct adapt_collector
{
  /** Collect adapt actions for all elements in the forest.
     * \param [in] forest_from        The forest containing the elements to be adapted.
     * \param [out] actions           The vector to store the adapt actions for all elements.
     * \param [in] callback           The callback function to determine the adapt action for each element.
     */
  void
  collect_actions (const t8_forest_t forest_from, std::vector<t8_adapt::action> &actions,
                   t8_adapt::element_callback callback)
  {
    T8_ASSERT (forest_from != nullptr);

    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
    actions.resize (static_cast<size_t> (local_num_elements));

    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_element_array_t *tree_elements_from = &tree_from->leaf_elements;
      const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (tree_elements_from);
      T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
      const t8_eclass_t tree_class = tree_from->eclass;
      const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

      for (t8_locidx_t el_considered = 0; el_considered < num_el_from; el_considered++) {
        const t8_element_t *element_from = t8_element_array_index_locidx (tree_elements_from, el_considered);
        actions[el_offset + el_considered]
          = callback (forest_from, ltree_id, el_considered, element_from, scheme, tree_class);
      }
      el_offset += num_el_from;
    }
  }
};

/**
 * Collect adaptation actions for all elements in the source forest.
 */
struct batched_adapt_collector
{
  /**
   * Collect adaptation actions for all elements in the source forest.
   * \param [in] forest_from      The source forest from which to collect adaptation actions.
   * \param [out] adapt_actions   The vector to store the collected adaptation actions.
   * \param [in] callback         The callback function to determine the adaptation action for each element.
   * 
   * \warning Currently only used for testing purposes.
   */
  void
  collect_actions (const t8_forest_t forest_from, std::vector<t8_adapt::action> &adapt_actions,
                   t8_adapt::batched_element_callback callback)
  {
    T8_ASSERT (forest_from != nullptr);

    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
    adapt_actions.resize (local_num_elements);

    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

    /* For each element get the adaptation action */
    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_eclass_t tree_class = tree_from->eclass;
      const t8_element_array_t *elements_from = &tree_from->leaf_elements;
      callback (forest_from, ltree_id, elements_from, scheme, tree_class, adapt_actions);
    }
  }
};

/** Standard family checker implementation. */
struct family_checker
{
  /**
     * Check if the elements in the given array are siblings0 in the tree and form a family.
     * \param [in] tree_elements_from The array of elements in the tree.
     * \param [in] offset              The offset to start checking from.
     * \param [in] scheme              The scheme to use for checking.
     * \param [in] tree_class          The class of the tree.
     * \return True if the elements are siblings and form a family, false otherwise.
     */
  bool
  family_check (const t8_element_array_t *tree_elements_from, const t8_locidx_t offset, const t8_scheme *scheme,
                const t8_eclass_t tree_class)
  {
    const t8_locidx_t num_elements_from = (t8_locidx_t) t8_element_array_get_count (tree_elements_from);
    const t8_element_t *first_element_from = t8_element_array_index_locidx (tree_elements_from, offset);
    const int num_siblings = scheme->element_get_num_siblings (tree_class, first_element_from);
    std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN> children;
    int added_siblings = 0;
    for (int isibling = 0; isibling < num_siblings && offset + (t8_locidx_t) isibling < num_elements_from; isibling++) {
      children[isibling] = const_cast<t8_element_t *> (
        t8_element_array_index_locidx (tree_elements_from, offset + (t8_locidx_t) isibling));
      if (scheme->element_get_child_id (tree_class, children[isibling]) != isibling) {
        return false;
      }
      added_siblings++;
    }
    if (added_siblings != num_siblings) {
      return false;
    }
    const bool is_family = scheme->elements_are_family (tree_class, children.data ());
    return is_family;
  };
};

/** Function to manipulate elements based on the specified adaptation action.
     * \tparam action The adaptation action to be performed.
     * \param [in, out] elements         The element array to be modified.
     * \param [in]     elements_from    The source element array.
     * \param [in]     scheme           The element scheme.
     * \param [in]     tree_class       The eclass of the tree used by the scheme
     * \param [in]     elements_index   The index in the target element array.
     * \param [in]     elements_from_index The index in the source element array.
     * \return                        The number of elements created in the target array.
     */
template <int action>
t8_locidx_t
manipulate_elements (t8_element_array_t *elements, [[maybe_unused]] const t8_element_array_t *const elements_from,
                     const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t elements_index,
                     const t8_locidx_t elements_from_index);


/** Specialization for the KEEP action. No element modification is performed; 
     * the element is copied as is.
     * \see manipulate_elements
     */
template <>
t8_locidx_t
manipulate_elements<t8_adapt::action::KEEP> (t8_element_array_t *elements,
                                             [[maybe_unused]] const t8_element_array_t *const elements_from,
                                             const t8_scheme *scheme, const t8_eclass_t tree_class,
                                             const t8_locidx_t elements_index,
                                             [[maybe_unused]] const t8_locidx_t elements_from_index)
{
  t8_element_t *element = t8_element_array_push (elements);
  const t8_element_t *element_from = t8_element_array_index_locidx (elements, elements_index);
  scheme->element_copy (tree_class, element_from, element);
  return 1;
};

/** Specialization for the COARSEN action. The parent element of the source elements is created 
     * in the target array.
     * \see manipulate_elements
     */
template <>
t8_locidx_t
manipulate_elements<t8_adapt::action::COARSEN> (t8_element_array_t *elements,
                                                [[maybe_unused]] const t8_element_array_t *const elements_from,
                                                const t8_scheme *scheme, const t8_eclass_t tree_class,
                                                const t8_locidx_t elements_index,
                                                [[maybe_unused]] const t8_locidx_t elements_from_index)
{
  t8_element_t *element = t8_element_array_push (elements);
  const t8_element_t *element_from = t8_element_array_index_locidx (elements, elements_index);
  T8_ASSERT (scheme->element_get_level (tree_class, element_from) > 0);
  scheme->element_get_parent (tree_class, element_from, element);

  /* Hier eventuell noch was mit num_children = num_siblings*/
  return 1;
};

/** Specialization for the REFINE action. The children elements of the source element are created 
     * in the target array.
     * \see manipulate_elements
     */
template <>
t8_locidx_t
manipulate_elements<t8_adapt::action::REFINE> (t8_element_array_t *elements,
                                               const t8_element_array_t *const elements_from, const t8_scheme *scheme,
                                               const t8_eclass_t tree_class, const t8_locidx_t elements_index,
                                               const t8_locidx_t elements_from_index)
{
  const t8_element_t *element_from = t8_element_array_index_locidx (elements_from, elements_from_index);
  const int num_children = scheme->element_get_num_children (tree_class, element_from);
  (void) t8_element_array_push_count (elements, num_children);
  std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN> children;
  for (int ichildren = 0; ichildren < num_children; ichildren++) {
    children[ichildren] = t8_element_array_index_locidx_mutable (elements, elements_index + ichildren);
  }
  scheme->element_get_children (tree_class, element_from, num_children, children.data ());
  return num_children;
};


t8_locidx_t
recursive_refine (t8_element_array_t *elements, t8_element_array_t *const elements_from,
                  const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t el_refine,
                  [[maybe_unused]]const t8_locidx_t el_offset, [[maybe_unused]]t8_locidx_t &el_inserted,
                  [[maybe_unused]]const std::vector<t8_adapt::action> &actions, [[maybe_unused]]const bool is_family, 
                  t8_adapt::element_callback callback, const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_locidx_t lelement_id)
{
  std::vector<t8_element_t *> element_stack;
  t8_element_t *current_element = t8_element_array_index_locidx_mutable (elements_from, el_refine);
  element_stack.push_back(current_element);
  t8_element_t **el_buffer;
  t8_locidx_t num_el_added = 0;

  while (!element_stack.empty()){
    el_buffer[0] = element_stack.back();
    element_stack.pop_back();
    const int num_children = scheme->element_get_num_children(tree_class, el_buffer[0]);
    /* TODO: Is lelement_id right in this case? */
    const t8_adapt::action refine = callback (forest, ltreeid, lelement_id, el_buffer[0], scheme, tree_class);
    T8_ASSERT (refine != t8_adapt::action::COARSEN);
    if (refine == t8_adapt::action::REFINE)
    {
      if (scheme->element_get_level(tree_class, el_buffer[0]) < forest->maxlevel &&
          scheme->element_is_refinable (tree_class, el_buffer[0])) {
            scheme->element_new (tree_class, num_children-1, el_buffer + 1);
            scheme->element_get_children(tree_class, el_buffer[0], num_children, el_buffer);
            element_stack.reserve (element_stack.size() + num_children);
            for (int ci = num_children - 1; ci >= 0; ci--)
            {
              element_stack.insert(element_stack.begin(), el_buffer[ci]);
            }
      }
    }
    else {
      T8_ASSERT(refine == t8_adapt::action::KEEP);
      num_el_added++;
      t8_element *el_insert = t8_element_array_push(elements);
      scheme->element_copy(tree_class, el_buffer[0], el_insert);
      scheme->element_destroy (tree_class, 1, el_buffer);
    }
  }
  return num_el_added;
};


t8_locidx_t
recursive_coarsen (t8_element_array_t *elements, [[maybe_unused]]const t8_element_array_t *const elements_from,
                  const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t el_coarsen,
                  [[maybe_unused]]const t8_locidx_t el_offset, t8_locidx_t &el_inserted,
                  [[maybe_unused]]const std::vector<t8_adapt::action> &actions, bool is_family, 
                  t8_adapt::element_callback callback, const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_locidx_t lelement_id)
{
  T8_ASSERT (el_coarsen >= 0);
  const t8_element_t *element = t8_element_array_index_locidx (elements, el_inserted-1);
  T8_ASSERT (scheme->element_get_level (tree_class, element) > 0);
  int num_siblings = scheme->element_get_num_siblings (tree_class, element);
  t8_locidx_t elements_in_array = t8_element_array_get_count(elements);
  t8_element_t **family = T8_ALLOC (t8_element_t *, num_siblings);
  t8_locidx_t position = t8_forest_pos (forest, tree_class, scheme, elements, el_inserted - 1);
  while (is_family && position >= el_coarsen && position < elements_in_array) {
    for (t8_locidx_t ielement = 0;ielement < (t8_locidx_t) num_siblings && position + ielement < elements_in_array; ielement++) {
      family[ielement] = t8_element_array_index_locidx_mutable(elements, position + ielement);
    }
    int num_elements_to_callback;
    if (scheme->elements_are_family (tree_class, family)) {
      num_elements_to_callback = num_siblings;
    }
    else {
      is_family = false;
      num_elements_to_callback = 1;
    }
    t8_adapt::action cb_result = t8_adapt::action::COARSEN;
    for (int ielement = 0; ielement < num_elements_to_callback; ielement++){

      cb_result = callback (forest, ltreeid, lelement_id, family[ielement], scheme, tree_class) != t8_adapt::action::COARSEN ? t8_adapt::action::REFINE : t8_adapt::action::COARSEN;
    }
    if (is_family &&
      cb_result == t8_adapt::action::COARSEN)
    {
      el_inserted -= (t8_locidx_t)(num_elements_to_callback - 1);
      scheme->element_get_parent(tree_class, family[0], family[0]);
      elements_in_array -= (t8_locidx_t)(num_elements_to_callback - 1);
      num_siblings = scheme->element_get_num_siblings(tree_class, family[0]);
      t8_element_array_resize(elements, elements_in_array);
      T8_ASSERT(el_inserted - 1 == position);
      position = t8_forest_pos(forest, tree_class, scheme, elements, position);
    }
    else {
      is_family = false;
    }
  }
  T8_FREE(family);
  return elements_in_array;
};
/** Standard element manipulator class for adapting elements between trees. */
struct manipulator
{
  void
  recursive_manipulate (t8_element_array_t *elements, const t8_element_array_t *const elements_from,
                       const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t el_coarsen,
                       const t8_locidx_t el_offset, t8_locidx_t &el_inserted,
                       const std::vector<t8_adapt::action> &actions, const bool is_family)
  {
    t8_adapt::action iaction = actions[el_coarsen];
    if (!is_family && iaction == t8_adapt::action::COARSEN) {
      iaction = t8_adapt::action::KEEP;
    }
    /* Check that all siblings want to be coarsened */
    if (is_family && iaction == t8_adapt::action::COARSEN) {
      const t8_element_t *current_element = t8_element_array_index_locidx (elements_from, el_coarsen);
      const int num_siblings = scheme->element_get_num_siblings (tree_class, current_element);
      const auto start = actions.begin () + static_cast<size_t> (el_coarsen);
      const auto end = start + static_cast<size_t> (num_siblings);
      if (!std::all_of (start, end, [] (const t8_adapt::action &action) { return action == t8_adapt::action::COARSEN; })) {
        iaction = t8_adapt::action::KEEP;
      }
    }
    switch (static_cast<int> (iaction)) {
    case t8_adapt::action::COARSEN:
      el_inserted += manipulate_elements<t8_adapt::action::COARSEN> (elements, elements_from, scheme, tree_class,
                                                                     el_offset + el_inserted, el_coarsen);
      break;
    case t8_adapt::action::KEEP:
      el_inserted += manipulate_elements<t8_adapt::action::KEEP> (elements, elements_from, scheme, tree_class,
                                                                  el_offset + el_inserted, el_coarsen);
      break;
    case t8_adapt::action::REFINE:
      el_inserted += manipulate_elements<t8_adapt::action::REFINE> (elements, elements_from, scheme, tree_class,
                                                                    el_offset + el_inserted, el_coarsen);
      break;
    default: {
      t8_errorf ("Unknown adapt action.\n");
      SC_ABORT_NOT_REACHED ();
      break;
    }
    }
  };
          
  /** Manipulate elements based on the given adapt action.
         * \param [in,out] elements          The array of elements to be manipulated.
         * \param [in] elements_from         The array of elements from the source tree.
         * \param [in] scheme                The scheme to use for manipulation.
         * \param [in] tree_class            The class of the tree.
         * \param [in] el_considered         The index of the element being considered.
         * \param [in] el_offset             The offset of the element being considered.
         * \param [in,out] el_inserted       The index of the elements in the target tree.
         * \param [in] actions               The global adapt actions vector.
         * \param [in] is_family             Whether the current elements form a family.
         */
  void
  element_manipulator (t8_element_array_t *elements, const t8_element_array_t *const elements_from,
                       const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t el_considered,
                       const t8_locidx_t el_offset, t8_locidx_t &el_inserted,
                       const std::vector<t8_adapt::action> &actions, const bool is_family)
  {
    t8_adapt::action iaction = actions[el_considered];
    if (!is_family && iaction == t8_adapt::action::COARSEN) {
      iaction = t8_adapt::action::KEEP;
    }
    /* Check that all siblings want to be coarsened */
    if (is_family && iaction == t8_adapt::action::COARSEN) {
      const t8_element_t *current_element = t8_element_array_index_locidx (elements_from, el_considered);
      const int num_siblings = scheme->element_get_num_siblings (tree_class, current_element);
      const auto start = actions.begin () + static_cast<size_t> (el_considered);
      const auto end = start + static_cast<size_t> (num_siblings);
      if (!std::all_of (start, end, [] (const t8_adapt::action &action) { return action == t8_adapt::action::COARSEN; })) {
        iaction = t8_adapt::action::KEEP;
      }
    }

    switch (static_cast<int> (iaction)) {
    case t8_adapt::action::COARSEN:
      el_inserted += manipulate_elements<t8_adapt::action::COARSEN> (elements, elements_from, scheme, tree_class,
                                                                     el_offset + el_inserted, el_considered);
      break;
    case t8_adapt::action::KEEP:
      el_inserted += manipulate_elements<t8_adapt::action::KEEP> (elements, elements_from, scheme, tree_class,
                                                                  el_offset + el_inserted, el_considered);
      break;
    case t8_adapt::action::REFINE:
      el_inserted += manipulate_elements<t8_adapt::action::REFINE> (elements, elements_from, scheme, tree_class,
                                                                    el_offset + el_inserted, el_considered);
      break;
    default: {
      t8_errorf ("Unknown adapt action.\n");
      SC_ABORT_NOT_REACHED ();
      break;
    }
    }
  };
};
};  // namespace t8_standard_adapt
