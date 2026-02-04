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

#ifndef T8_FOREST_ADAPT_BATCHED_HXX
#define T8_FOREST_ADAPT_BATCHED_HXX

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>

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
template <t8_forest_adapt_namespace::adapt_action action>
t8_locidx_t
manipulate_elements (t8_element_array_t *elements, const t8_element_array_t *const elements_from,
                     const t8_scheme *scheme, const t8_eclass_t tree_class, const t8_locidx_t elements_index,
                     const t8_locidx_t elements_from_index);

/** Specialization for the KEEP action. No element modification is performed; 
   * the element is copied as is.
   * \see manipulate_elements
  */
template <>
t8_locidx_t
manipulate_elements<t8_forest_adapt_namespace::adapt_action::KEEP> (
  t8_element_array_t *elements, const t8_element_array_t *const elements_from, const t8_scheme *scheme,
  const t8_eclass_t tree_class, const t8_locidx_t elements_index, const t8_locidx_t elements_from_index)
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
manipulate_elements<t8_forest_adapt_namespace::adapt_action::COARSEN> (
  t8_element_array_t *elements, const t8_element_array_t *const elements_from, const t8_scheme *scheme,
  const t8_eclass_t tree_class, const t8_locidx_t elements_index, const t8_locidx_t elements_from_index)
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
manipulate_elements<adapt_action::REFINE> (t8_element_array_t *elements, const t8_element_array_t *const elements_from,
                                           const t8_scheme *scheme, const t8_eclass_t tree_class,
                                           const t8_locidx_t elements_index, const t8_locidx_t elements_from_index)
{
  const t8_element_t *element_from = t8_element_array_index_locidx (elements_from, elements_from_index);
  const int num_children = scheme->element_get_num_children (tree_class, element_from);
  /* CONTINUE WORK HERE */
  (void) t8_element_array_push_count (elements, num_children);
  std::vector<t8_element_t *> children (num_children);
  for (int ichildren = 0; ichildren < num_children; ichildren++) {
    children[ichildren] = t8_element_array_index_locidx_mutable (elements, elements_index + ichildren);
  }
  scheme->element_get_children (tree_class, element_from, num_children, children.data ());
  return num_children;
};

/**
     * Collect adaptation actions for all elements in the source forest.
     */
struct batched_adapt_collector
{
  void
  collect_adapt_actions (const t8_forest_t forest_from, std::vector<adapt_action> &adapt_actions,
                         element_callback callback)
  {
    T8_ASSERT (forest_from != nullptr);

    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
    adapt_actions.resize (local_num_elements);

    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

    /* For each element get the adaptation action */
    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_eclass_t tree_class = tree_from->eclass;
      const t8_element_array_t *elements_from = &tree_from->leaf_elements;
      const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (elements_from);
      for (t8_locidx_t i = 0; i < num_el_from; i++) {
        const t8_element_t *element = t8_element_array_index_locidx (elements_from, i);
        adapt_actions[el_offset + i] = callback (forest_from, ltree_id, element, scheme, tree_class);
      }
      el_offset += num_el_from;
    }
  }
};

#endif /* T8_FOREST_ADAPT_BATCHED_HXX */
