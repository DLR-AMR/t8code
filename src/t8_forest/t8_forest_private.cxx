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

#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_general.h>
#include <algorithm>

T8_EXTERN_C_BEGIN ();

const t8_element_t *
t8_forest_get_tree_leaf_element (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= elem_in_tree && elem_in_tree < t8_forest_get_tree_leaf_element_count (tree));
  return t8_element_array_index_locidx (&tree->leaf_elements, elem_in_tree);
}

t8_element_t *
t8_forest_get_tree_leaf_element_mutable (t8_tree_t tree, t8_locidx_t elem_in_tree)
{
  return (t8_element_t *) t8_forest_get_tree_leaf_element (tree, elem_in_tree);
}

const t8_element_array_t *
t8_forest_get_tree_leaf_element_array (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  return &t8_forest_get_tree (forest, ltreeid)->leaf_elements;
}

t8_element_array_t *
t8_forest_get_tree_leaf_element_array_mutable (const t8_forest_t forest, t8_locidx_t ltreeid)
{
  return (t8_element_array_t *) t8_forest_get_tree_leaf_element_array (forest, ltreeid);
}

T8_EXTERN_C_END ();

/* TODO: does the search fail when element_level is smaller then levels in the array?
         For example entering the search with the root element or a level 1 element
         and the array contains much finer elements.
         Will it still return the largest index, or just any index?
 */
/* TODO: This may be implementable with std::partition_point, which would yield an easier implementation.
         Need to check.
 */
/** \brief Search for a linear element id (at level element_level) in a sorted array of
 * elements. If the element does not exist, return the largest index i
 * such that the element at position i has a smaller id than the given one.
 * If no such i exists, return -1.
 */
t8_locidx_t
t8_forest_bin_search_lower (const t8_element_array_view &elements, const t8_linearidx_t element_id,
                            const int element_level)
{
  const t8_scheme *scheme = elements.get_scheme ();
  const t8_eclass_t tree_class = elements.get_tree_class ();

  /* We search for the first element in the array that is greater than the given element id. */
  const t8_element_array_view::iterator elem_iter = std::upper_bound (
    elements.begin (), elements.end (), element_id,
    [element_level, scheme, tree_class] (const t8_linearidx_t element_id_, const t8_element_t *element) {
      return element_id_ < scheme->element_get_linear_id (tree_class, element, element_level);
    });

  /* All elements before the found one have a smaller or equal id, so we return the index before it.
   * If the found element is the first one (or the array is empty), no element has a smaller or equal id
   * and we return -1. If no element has a greater id, the search returns the end iterator and we return
   * the last index of the array. */
  return static_cast<t8_locidx_t> (elem_iter - elements.begin ()) - 1;
}

/** \brief Search for a linear element id (at level element_level) in a sorted array of
 * elements. If the element does not exist, return the smallest index i
 * such that the element at position i has a larger id than the given one.
 * If no such i exists, return -1.
 */
t8_locidx_t
t8_forest_bin_search_upper (const t8_element_array_view &elements, const t8_linearidx_t element_id,
                            const int element_level)
{
  const t8_scheme *scheme = elements.get_scheme ();
  const t8_eclass_t tree_class = elements.get_tree_class ();

  /* We search for the first element E in the array, where ID(E) < element_id is false.
     Thus, E is the first element with ID(E) >= element_id . */
  const t8_element_array_view::iterator elem_iter = std::lower_bound (
    elements.begin (), elements.end (), element_id,
    [element_level, scheme, tree_class] (const t8_element_t *element, const t8_linearidx_t element_id_) {
      return scheme->element_get_linear_id (tree_class, element, element_level) < element_id_;
    });

  /* In case no element has a greater or equal id (or the array is empty), the binary search returns
   * the end iterator. */
  if (elem_iter == elements.end ()) {
    // No element was found.
    return -1;
  }
  return static_cast<t8_locidx_t> (elem_iter - elements.begin ());
}

/** \brief Search for a linear element id (at level element_level) in a sorted array of
 * elements. If the element does not exist, return the first index i such that
 * the element at position i is an ancestor or descendant of the element corresponding to the element id.
 * If no such i exists, return -1.
 */
t8_locidx_t
t8_forest_bin_search_first_descendant_ancestor (const t8_element_array_view &elements, const t8_element_t *element,
                                                const t8_element_t **element_found)
{
  /* This search works as follows:

  Let E denote the element with element_id at level L.
  If an ancestor or descendant of E exists in the array then they are either:
    A: The search result of t8_forest_bin_search_lower
    B: The search result of t8_forest_bin_search_upper

  Let ID(element,level) denote the linear id of an element at a given level.

  Case A: There is an element F in the array that is E itself or an ancestor of E (i.e. level(F) <= level(E)).
          In that case
          ID(E,L) >= ID(F,L) and there can be no element with id in between (since it would also be an ancestor of E).
          Then F will be the search result of t8_forest_bin_search_lower
  Case B: There is an element F in the array that is a descendant of E and it has the smallest index in the array of all descendants.
          Then
          ID(E,L) = ID(F,L)
          and also
          ID(E,L) = ID(D,L) for all other descendants of E.
          But since F is the first it will be the search result of t8_forest_bin_search_upper.
  Case C: There is no descendant or ancestor of E in the array. In both cases t8_forest_bin_search_lower and
          t8_forest_bin_search_upper may find elements but the results will not be ancestors/descendants of E.

   From this, we determine the following algorithm:

    1. Query t8_forest_bin_search_lower with N.
    2. If no element was found, or the resulting element is not an ancestor of N.
    3. Query t8_forest_bin_search_upper with N.
    3. If an element was found and it is a descendant of N, we found our element.
    4. If not, no element was found.
  */

  /* Compute the element's level and linear id. In order to do so,
   * we first need the scheme and eclass. */
  const t8_scheme *scheme = elements.get_scheme ();
  const t8_eclass_t eclass = elements.get_tree_class ();
  const int element_level = scheme->element_get_level (eclass, element);
  const t8_linearidx_t element_id = scheme->element_get_linear_id (eclass, element, element_level);

  const t8_locidx_t search_pos_lower = t8_forest_bin_search_lower (elements, element_id, element_level);

  /* Get the element at the current position. */
  if (search_pos_lower >= 0) {
    *element_found = elements[search_pos_lower];
    const bool is_ancestor = scheme->element_is_ancestor (eclass, *element_found, element);
    if (is_ancestor) {
      /* The element at this position is an ancestor or descendant. */
      return search_pos_lower;
    }
  }
  /* t8_forest_bin_search_lower did not return a result or an ancestor. */

  const t8_locidx_t search_pos_upper = t8_forest_bin_search_upper (elements, element_id, element_level);
  if (search_pos_upper >= 0) {
    *element_found = elements[search_pos_upper];
    const bool is_descendant = scheme->element_is_ancestor (eclass, element, *element_found);
    if (is_descendant) {
      /* The element at this position is an ancestor or descendant. */
      return search_pos_upper;
    }
  }
  // No ancestor or descendant was found
  *element_found = nullptr;
  return -1;
}
