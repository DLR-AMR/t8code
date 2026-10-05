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

/** \file t8_gtest_adapt_callbacks.cxx
* Provide forest adapt callback functions that we use in our tests.
*/

#include <test/t8_gtest_adapt_callbacks.hxx>

/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest.
 * 
 * This adapt callbacks requires an integer as forest user data.
 * This integer is the maximum refinement level.
 */
int
t8_test_adapt_first_child (const t8_forest_t forest, [[maybe_unused]]const t8_locidx_t which_tree, [[maybe_unused]]const t8_locidx_t lelement_id,
                         const t8_element_t *element, const t8_scheme_c *scheme, const t8_eclass_t tree_class)
{
  const int level = scheme->element_get_level (tree_class, element);

  /* we set a maximum refinement level as forest user data */
  const int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  const int child_id = scheme->element_get_child_id (tree_class, element);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}
