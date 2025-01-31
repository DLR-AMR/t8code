/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>

/*
 * In this file we test whether the t8_element_count_leaves{_from_root}
 * function return the correct values for the default scheme.
 * This value should be  2^(dim * (level - element_level)) for
 * level >= element_level (for eclass != T8_ECLASS_PYRAMID).
 * For level < element_level the value should be zero.
 */

/* Tests whether the leaf count for one additional level matches the number of children */

class class_element_leaves: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    eclass = scheme->get_eclass_scheme_eclass (eclass);
  }
  void
  TearDown () override
  {
    scheme->unref ();
  }
  const t8_scheme *scheme;
  t8_eclass_t eclass;
};

TEST_P (class_element_leaves, test_element_count_leaves_root)
{
  const int maxlevel = scheme->get_maxlevel (eclass);
  t8_gloidx_t compare_value = 1;
  t8_gloidx_t sum1 = 1;
  t8_gloidx_t sum2 = 1;

  for (int level = 0; level <= maxlevel; ++level) {
    const t8_gloidx_t leaf_count = scheme->count_leaves_from_root (eclass, level);
    ASSERT_EQ (leaf_count, compare_value)
      << "Incorrect leaf count " << leaf_count << " at eclass " << t8_eclass_to_string[eclass] << " and level " << level
      << " (expecting " << compare_value << ")";
    /* Multiply the compare_value with 2^dim (= number of children per element) */
    if (eclass == T8_ECLASS_PYRAMID) {
      sum1 *= 8;
      sum2 *= 6;
      compare_value = 2 * sum1 - sum2;
    }
    else {
      compare_value *= 1 << t8_eclass_to_dimension[eclass];
    }
  }
}

/* Tests whether the leaf count for the same level is equal to 1
 * and for smaller levels is 0 */
TEST_P (class_element_leaves, test_element_count_leaves_less_level)
{
  t8_element_t *element;
  const int maxlevel = scheme->get_maxlevel (eclass);

  /* Allocate memory for an element */
  scheme->element_new (eclass, 1, &element);
  for (int level = 0; level <= maxlevel; ++level) {
    /* Create the first element on this level */
    scheme->element_set_linear_id (eclass, element, level, 0);
    /* Count the leaves of this element */
    const t8_gloidx_t leaf_count_same_level = scheme->element_count_leaves (eclass, element, level);
    /* Check if equals 1 */
    ASSERT_EQ (leaf_count_same_level, 1);
    int lower_levels;
    for (lower_levels = level - 1; lower_levels >= 0; --lower_levels) {
      /* Count the leaves of this element on the lower levels */
      const t8_gloidx_t leaf_count = scheme->element_count_leaves (eclass, element, lower_levels);
      /* Check if equals 0 */
      ASSERT_EQ (leaf_count, 0) << "Incorrect leaf count " << leaf_count << " at eclass " << t8_eclass_to_string[eclass]
                                << " and level " << level << " for element level " << lower_levels << "(expecting 0)";
    }
  }
  /* Free the element's memory */
  scheme->element_destroy (eclass, 1, &element);
}

/* Tests whether the leaf count for one additional level matches the number of children */
TEST_P (class_element_leaves, test_element_count_leaves_one_level)
{
  t8_element_t *element;
  const int maxlevel = scheme->get_maxlevel (eclass);

  scheme->element_new (eclass, 1, &element);
  for (int level = 1; level < maxlevel; ++level) {
    /* Create the first element on the previous level */
    scheme->element_set_linear_id (eclass, element, level - 1, 0);
    /* Count the leaves of this element */
    const t8_gloidx_t leaf_count = scheme->element_count_leaves (eclass, element, level);
    /* Compute the number of children of the element */
    const int number_of_children = scheme->element_get_num_children (eclass, element);
    /* Check both values for equality */
    ASSERT_EQ (leaf_count, number_of_children)
      << "Incorrect leaf count " << leaf_count << " at eclass " << t8_eclass_to_string[eclass] << " and level " << level
      << " (expecting " << number_of_children << ")";
  }
  scheme->element_destroy (eclass, 1, &element);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_count_leaves, class_element_leaves, AllSchemes, print_all_schemes);
