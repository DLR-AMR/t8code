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

/** \file t8_gtest_adapt_callbacks.hxx
* Provide forest adapt callback functions that we use in our tests.
*/

#ifndef T8_GTEST_ADAPT_CALLBACKS
#define T8_GTEST_ADAPT_CALLBACKS

#include <src/t8_forest/t8_forest_general.h>

/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest.
 * 
 * This adapt callbacks requires an integer as forest user data.
 * This integer is the maximum refinement level.
 */
int
t8_test_adapt_first_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_eclass_scheme_c *scheme, const int is_family, const int num_elements,
                           t8_element_t *elements[]);

#endif /* T8_GTEST_ADAPT_CALLBACKS */
