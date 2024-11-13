/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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

#include "t8_forest/t8_forest_search/t8_forest_search.hxx"
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element.hxx>

void
t8_search_base::search_recursion (const t8_locidx_t ltreeid, t8_element_t *element, const t8_eclass_scheme_c *ts,
                                  t8_element_array_t *leaf_elements, const t8_locidx_t tree_lindex_of_first_leaf)
{
  /* Assertions to check for necessary requirements */
  /* The forest must be committed */
  T8_ASSERT (t8_forest_is_committed (this->forest));
  /* The tree must be local */
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (this->forest));

  const size_t elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }

  if (this->stop_due_to_queries ()) {
    return;
  }

  bool is_leaf = false;
  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);

    SC_CHECK_ABORT (ts->t8_element_level (element) <= ts->t8_element_level (leaf),
                    "Search: element level greater than leaf level\n");
    if (ts->t8_element_level (element) == ts->t8_element_level (leaf)) {
      T8_ASSERT (t8_forest_element_is_leaf (this->forest, leaf, ltreeid));
      T8_ASSERT (ts->t8_element_equal (element, leaf));
      /* The element is the leaf */
      is_leaf = true;
    }
  }
  /* Call the callback function for the element */
  const bool ret = check_element (ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf);

  if (!ret) {
    /* The function returned false. We abort the recursion */
    return;
  }
  std::vector<size_t> new_active_queries;
  this->check_queries (new_active_queries, ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf);

  if (is_leaf) {
    return;
  }

  /* Enter the recursion (the element is definitely not a leaf at this point) */
  /* We compute all children of E, compute their leaf arrays and call search_recursion */
  /* allocate the memory to store the children */
  const int num_children = ts->t8_element_num_children (element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  ts->t8_element_new (num_children, children);
  /* Memory for the indices that split the leaf_elements array */
  size_t *split_offsets = T8_ALLOC (size_t, num_children + 1);
  /* Compute the children */
  ts->t8_element_children (element, num_children, children);
  /* Split the leaves array in portions belonging to the children of element */
  t8_forest_split_array (element, leaf_elements, split_offsets);
  for (int ichild = 0; ichild < num_children; ichild++) {
    /* Check if there are any leaf elements for this child */
    const size_t indexa = split_offsets[ichild];     /* first leaf of this child */
    const size_t indexb = split_offsets[ichild + 1]; /* first leaf of next child */
    if (indexa < indexb) {
      t8_element_array_t child_leaves;
      /* There exist leaves of this child in leaf_elements,
       * we construct an array of these leaves */
      t8_element_array_init_view (&child_leaves, leaf_elements, indexa, indexb - indexa);
      /* Enter the recursion */
      search_recursion (ltreeid, children[ichild], ts, &child_leaves, indexa + tree_lindex_of_first_leaf);
      update_queries (new_active_queries);
    }
  }

  /* clean-up */
  ts->t8_element_destroy (num_children, children);
  T8_FREE (children);
  T8_FREE (split_offsets);
}

void
t8_search_base::search_tree (const t8_locidx_t ltreeid)
{
  const t8_eclass_t eclass = t8_forest_get_eclass (this->forest, ltreeid);
  const t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (this->forest, eclass);
  t8_element_array_t *leaf_elements = t8_forest_tree_get_leaves (this->forest, ltreeid);

  /* assert for empty tree */
  T8_ASSERT (t8_element_array_get_count (leaf_elements) >= 0);
  /* Get the first and last leaf of this tree */
  const t8_element_t *first_el = t8_element_array_index_locidx (leaf_elements, 0);
  const t8_element_t *last_el
    = t8_element_array_index_locidx (leaf_elements, t8_element_array_get_count (leaf_elements) - 1);
  /* Compute their nearest common ancestor */
  t8_element_t *nca;
  ts->t8_element_new (1, &nca);
  ts->t8_element_nca (first_el, last_el, nca);

  /* Start the top-down search */
  this->search_recursion (ltreeid, nca, ts, leaf_elements, 0);

  ts->t8_element_destroy (1, &nca);
}

void
t8_search_base::do_search ()
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (this->forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    this->search_tree (itree);
  }
}
