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

/**
 * \file This file defines some helper functions used for the partition-for-coarsening feature.
*/

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_schemes/t8_scheme.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_element.h>
#include <algorithm>

t8_locidx_t
t8_forest_pfc_extreme_local_sibling (const t8_scheme_c *scheme, const t8_tree_t tree,
                                     const t8_locidx_t start_element_id_in_tree, const int signed_increment)
{
  T8_ASSERT (signed_increment == 1 || signed_increment == -1);

  // Initialization and memory allocation.
  const t8_eclass_t tree_class = tree->eclass;
  t8_element_t *parent_possible_sibling, *parent_start;
  t8_element_new (scheme, tree_class, 1, &parent_possible_sibling);
  t8_element_new (scheme, tree_class, 1, &parent_start);

  // Determine start element from tree and start ID within tree.
  const t8_element_t *start_element = t8_forest_get_tree_leaf_element (tree, start_element_id_in_tree);

  // If the start element is of level zero, i.e., the root, it does not have any siblings.
  if (scheme->element_get_level (tree_class, start_element) == 0) {
    return start_element_id_in_tree;
  }

  // Get parent of start element.
  scheme->element_get_parent (tree_class, start_element, parent_start);

  // Determine the parent's number of children.
  const int num_children = scheme->element_get_num_children (tree_class, parent_start);

  // Determine increment and bound to be used within element loop:
  // (a) increment = -1 and lower bound, or
  // (b) increment = +1 and upper bound.
  t8_locidx_t extreme_check_id_in_tree;
  if (signed_increment < 0) {
    extreme_check_id_in_tree = SC_MAX (0, start_element_id_in_tree - num_children);
  }
  else {
    extreme_check_id_in_tree
      = SC_MIN (start_element_id_in_tree + num_children, t8_forest_get_tree_leaf_element_count (tree)) - 1;
  }

  // Initialize extreme_sibling_id_in_tree.
  t8_locidx_t extreme_sibling_id_in_tree = start_element_id_in_tree;

  // Loop over local IDs of all elements that may form a family
  for (t8_locidx_t ielem = start_element_id_in_tree; (ielem - extreme_check_id_in_tree) * signed_increment <= 0;
       ielem += signed_increment) {

    // Get element from iteration index.
    const t8_element_t *possible_sibling = t8_forest_get_tree_leaf_element (tree, ielem);

    // Determine parent and check whether it matches parent_start:
    // - if it does, the current element is a sibling of start_element. Thus, extreme_sibling_id_in_tree is updated.
    scheme->element_get_parent (tree_class, possible_sibling, parent_possible_sibling);
    if (scheme->element_is_equal (tree_class, parent_start, parent_possible_sibling)) {
      extreme_sibling_id_in_tree = ielem;
    }
    // - else, the iteration has left the family and we can exit.
    else {
      break;
    }
  }

  // Deallocation
  t8_element_destroy (scheme, tree_class, 1, &parent_possible_sibling);
  t8_element_destroy (scheme, tree_class, 1, &parent_start);

  // Return extreme sibling ID.
  return extreme_sibling_id_in_tree;
}

void
t8_forest_pfc_helper_index_in_tree_from_globalid (const t8_forest_t forest, const t8_gloidx_t gelement_id,
                                                  t8_gloidx_t &gtree_id, t8_tree_t &tree, t8_locidx_t &index_in_tree,
                                                  t8_element_t *&element)
{
  // Determine local element ID by subtracting given and the process' first global ID.
  const t8_gloidx_t global_id_of_first_local_element = t8_forest_get_first_local_leaf_element_id (forest);
  const t8_locidx_t lelement_id = (t8_locidx_t) (gelement_id - global_id_of_first_local_element);

  // Determine the element (as pointer) and the local tree ID.
  t8_locidx_t ltree_id;
  element = t8_forest_get_leaf_element (forest, lelement_id, &ltree_id);

  // From local tree ID, get a pointer to the tree.
  tree = t8_forest_get_tree (forest, ltree_id);

  // Determine the index within the tree - and run sanity check.
  index_in_tree = lelement_id - tree->elements_offset;
  T8_ASSERT (element == t8_forest_get_tree_leaf_element (tree, index_in_tree));

  // Compute global tree ID as the local one plus the process' first tree ID.
  const t8_locidx_t first_local_tree_id = t8_forest_get_first_local_tree_id (forest);
  gtree_id = first_local_tree_id + ltree_id;
}
