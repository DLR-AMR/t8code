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
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/*
 * In this file we test whether the t8_element_count_leaves{_from_root}
 * function return the correct values for the default scheme.
 * This value should be  2^(dim * (level - element_level)) for
 * level >= element_level (for eclass != T8_ECLASS_PYRAMID).
 * For level < element_level the value should be zero.
 */

/* Tests whether the leaf count for one additional level matches the number of children */

struct class_element_leaves: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    eclass = scheme->get_eclass_scheme_eclass (eclass);

    scheme->element_new (eclass, 1, &element);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &element);
    scheme->unref ();
  }
  t8_element_t *element;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
};

static int
adapt_all ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
           [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
           [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme_c *scheme,
           [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t *elements[])
{
  return 1;
}

TEST_P (class_element_leaves, test_element_count_leaves_root)
{
#if T8CODE_TEST_LEVEL >= 1
  const int maxlevel = 4;
#else
  const int maxlevel = 6;
#endif
  t8_gloidx_t compare_value = 1;
  t8_gloidx_t test_value = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  t8_cmesh_ref (cmesh);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, sc_MPI_COMM_WORLD);
  scheme->ref ();
  for (int level = 0; level <= maxlevel; ++level) {
    const t8_gloidx_t leaf_count = scheme->count_leaves_from_root (eclass, level);
    ASSERT_EQ (leaf_count, compare_value)
      << "Incorrect leaf count " << leaf_count << " at eclass " << t8_eclass_to_string[eclass] << " and level " << level
      << " (expecting " << compare_value << ")";
    forest = t8_forest_new_adapt (forest, adapt_all, 0, 0, NULL);
    compare_value = t8_forest_get_local_num_leaf_elements (forest);

    if (!scheme->refines_irregular (eclass)) {
      test_value *= 1 << t8_eclass_to_dimension[eclass];
      EXPECT_EQ (test_value, compare_value);
    }
  }
  t8_cmesh_unref (&cmesh);
  t8_forest_unref (&forest);
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
