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
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_adapt_callbacks.hxx>
#include <test/t8_gtest_macros.hxx>

/* In this test we check the t8_forest_element_is_leaf function.
 * Iterating over all cmesh test cases, we creat a uniform and an adaptive forest.
 * For each forest, we check that for each leaf element t8_forest_element_is_leaf returns true
 * and that it returns false for the parent and the first child.
 */

/* Maximum uniform level for forest. */

#if T8_TEST_LEVEL_INT >= 1
#define T8_IS_LEAF_MAX_LVL 3
#else
#define T8_IS_LEAF_MAX_LVL 4
#endif

struct element_is_leaf: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int>>
{
 protected:
  void
  SetUp () override
  {
    /* Construct a cmesh */
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    const t8_eclass_t tree_class = std::get<1> (std::get<0> (GetParam ()));
    const int level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = t8_cmesh_new_from_class (tree_class, sc_MPI_COMM_WORLD);

    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    //const int maxlevel = t8_forest_get_maxlevel (forest);
    int maxlevel = 7;
    const int recursive_adapt = 1;
    forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, 0, &maxlevel);
  }

  void
  TearDown () override
  {
    if (forest != NULL) {
      t8_forest_unref (&forest);
    }
    if (forest_adapt != NULL) {
      t8_forest_unref (&forest_adapt);
    }
  }

  t8_forest_t forest { NULL };
  t8_forest_t forest_adapt { NULL };
  const t8_scheme *scheme;
};

struct element_is_leaf_hybrid: public testing::TestWithParam<int>
{
 protected:
  void
  SetUp () override
  {
    /* Construct a cmesh */
    const int scheme_id = GetParam ();
    scheme = create_from_scheme_id (scheme_id);
    t8_cmesh_t cmesh = t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);
    const int level = 0;
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    int maxlevel = 7;
    const int recursive_adapt = 1;
    forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, 0, &maxlevel);
  }

  void
  TearDown () override
  {
    if (forest != NULL) {
      t8_forest_unref (&forest);
    }
    if (forest_adapt != NULL) {
      t8_forest_unref (&forest_adapt);
    }
  }

  t8_forest_t forest { NULL };
  t8_forest_t forest_adapt { NULL };
  const t8_scheme *scheme;
};

static void
t8_test_element_is_leaf_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Allocate memory to build a non-leaf element. */
    t8_element_t *not_leaf;
    scheme->element_new (tree_class, 1, &not_leaf);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      EXPECT_TRUE (t8_forest_element_is_leaf (forest, leaf_element, itree));
      /* Compute parent and first child of element and check that they are not in the tree */
      const int element_level = scheme->element_get_level (tree_class, leaf_element);
      if (element_level > 0) {
        scheme->element_get_parent (tree_class, leaf_element, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
      }
      if (element_level < scheme->get_maxlevel (tree_class)) {
        scheme->element_get_child (tree_class, leaf_element, 0, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
      }
    }
    scheme->element_destroy (tree_class, 1, &not_leaf);
  }
}

TEST_P (element_is_leaf, element_is_leaf)
{
  t8_test_element_is_leaf_for_forest (forest);
}

TEST_P (element_is_leaf, element_is_leaf_adapt)
{
  t8_test_element_is_leaf_for_forest (forest_adapt);
}

TEST_P (element_is_leaf_hybrid, element_is_leaf)
{
  t8_test_element_is_leaf_for_forest (forest);
}

TEST_P (element_is_leaf_hybrid, element_is_leaf_adapt)
{
  t8_test_element_is_leaf_for_forest (forest_adapt);
}

/* Define a lambda to beatify gtest output for tuples <level, cmesh>.
 * This will set the correct level and cmesh name as part of the test case name. */
auto pretty_print_eclass_scheme_and_level
  = [] (const testing::TestParamInfo<std::tuple<std::tuple<int, t8_eclass_t>, int>> &info) {
      std::string scheme = t8_scheme_to_string[std::get<0> (std::get<0> (info.param))];
      std::string eclass = t8_eclass_to_string[std::get<1> (std::get<0> (info.param))];
      std::string level = std::string ("_level_") + std::to_string (std::get<1> (info.param));
      return scheme + "_" + eclass + level;
    };

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_leaf, element_is_leaf,
                          testing::Combine (AllSchemes, testing::Range (0, T8_IS_LEAF_MAX_LVL)),
                          pretty_print_eclass_scheme_and_level);

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_leaf_hybrid, element_is_leaf_hybrid, AllSchemeCollections, print_scheme);
