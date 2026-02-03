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
#include <t8_forest/t8_forest_private.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_adapt_callbacks.hxx>
#include <test/t8_gtest_macros.hxx>

/* In these tests we check the t8_forest_bin_search_lower, t8_forest_bin_search_upper
 * and t8_forest_bin_search_first_descendant_ancestor functions.
 * Iterating over all cmesh test cases, we create a uniform and an adaptive forest.
 * For each forest, we 
 * TODO: Fill out comment
 */

/* Maximum uniform level for forest. */

#if T8_TEST_LEVEL_INT >= 1
#define T8_IS_LEAF_MAX_LVL 3
#else
#define T8_IS_LEAF_MAX_LVL 4
#endif

// TODO: ADd this class to common headers since it is reused
class t8_bin_search_tester: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int>> {
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

    // Construct a uniform forest
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    int maxlevel = 7;
    const int recursive_adapt = 1;
    // Construct an adaptive forest
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

/** Test the t8_forest_bin_search_upper function.
 * We iterate through all elements of a forest.
 * For each element E of level L, we call t8_forest_bin_search_lower on the forest's element array and expect
 * the element to be found.
 * For each element we then build one case where we search for an element that is not contained but will
 * match a different element in the element array:
 *  We compute the linear Id of E at level L-1 (if it is >=0 ) and search for the L-1 element with this id.
 *  We expect E to be found since it fulfills that its id is >= than the id we search for.
 * We then build one case per tree where we search for an element that is not contained and will not match
 * any other element:
 *  We compute the linear Id of the last element in the tree and add 1 to it (if it is not zero).
 *  We expect the search to not find anything.
*/
static void
t8_test_forest_bin_search_upper (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_element_array_t *leafs = t8_forest_tree_get_leaf_elements (forest, itree);

    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      const int element_level = scheme->element_get_level (tree_class, leaf_element);
      const t8_linearidx_t element_id = scheme->element_get_linear_id (tree_class, leaf_element, element_level);

      /* Search for a linear element id in a sorted array of
      * elements. If the element does not exist, return the largest index i
      * such that the element at position i has a smaller id than the given one.
      * If no such i exists, return -1. */
      const t8_locidx_t search_index = t8_forest_bin_search_upper (leafs, element_id, element_level);
      // We expect the leaf element to be found at position ielement
      EXPECT_EQ (search_index, ielement) << "Found wrong position of leaf element. Expected: " << ielement
                                         << " got: " << search_index;

      // If we increase the level, we expect the element to not be found, but the search
      // should return the index of the original element.
      if (element_level < scheme->get_maxlevel (tree_class)) {
        t8_debugf ("Computing element for level %i, Max is %i\n", element_level, T8_DLINE_MAXLEVEL);
        // TODO: The maxlevel eventually should be element dependent. I know that an element dependent maxlevel function was developed in a different branch (by Sandro?)
        const t8_linearidx_t element_id_at_next_level
          = scheme->element_get_linear_id (tree_class, leaf_element, element_level + 1);

        const t8_locidx_t search_index
          = t8_forest_bin_search_upper (leafs, element_id_at_next_level, element_level + 1);
        // We expect the leaf element to be found at position ielement
        EXPECT_EQ (search_index, ielement)
          << "Found wrong position of level " << element_level + 1 << " leaf element with id "
          << element_id_at_next_level << ". Expected: " << ielement << " got: " << search_index;
      }

      // Construct an element that is definitely not in the array and
      // does not have an element of smaller id in the array. We expect -1 as return.
      // We take the first element of the forest and subtract 1 from its id.
      if (ielement == num_elements_in_tree - 1) {
        const t8_linearidx_t element_not_found_id = element_id + 1;
        // Double check for possible conversion error, if element_id == MAX_POSSIBLE_VALUE. In that case, the
        // test logic fails. Should we ever run into this case, we need to rewrite this test accordingly.
        SC_CHECK_ABORTF (element_not_found_id > 0, "Invalid element id %li\n", element_not_found_id);
        const t8_locidx_t search_index = t8_forest_bin_search_upper (leafs, element_not_found_id, element_level);
        EXPECT_EQ (search_index, -1) << "Wrong return value for element that should not be in array. Expectec -1.";
      }
    }
  }
}

/** Test the t8_forest_bin_search_lower function.
 * We iterate through all elements of a forest.
 * For each element E of level L, we call t8_forest_bin_search_lower on the forest's element array and expect
 * the element to be found.
 * For each element we then build one case where we search for an element that is not contained but will
 * match a different element in the element array:
 *  We compute the linear Id of E at level L+1 (if not exceeding the maxlevel) and search for the L+1 element with this id.
 *  We expect E to be found since it fulfills that its id is <= than the id we search for.
 * We then build one case per tree where we search for an element that is not contained and will not match
 * any other element:
 *  We compute the linear Id of the first element in the tree and subtract 1 from it (if it is not zero).
 *  We expect the search to not find anything.
*/
static void
t8_test_forest_bin_search_lower (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_element_array_t *leafs = t8_forest_tree_get_leaf_elements (forest, itree);

    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      const int element_level = scheme->element_get_level (tree_class, leaf_element);
      const t8_linearidx_t element_id = scheme->element_get_linear_id (tree_class, leaf_element, element_level);

      /* Search for a linear element id in a sorted array of
      * elements. If the element does not exist, return the largest index i
      * such that the element at position i has a smaller id than the given one.
      * If no such i exists, return -1. */
      const t8_locidx_t search_index = t8_forest_bin_search_lower (leafs, element_id, element_level);
      // We expect the leaf element to be found at position ielement
      EXPECT_EQ (search_index, ielement) << "Found wrong position of leaf element. Expected: " << ielement
                                         << " got: " << search_index;

      // If we increase the level, we expect the element to not be found, but the search
      // should return the index of the original element.
      if (element_level < scheme->get_maxlevel (tree_class)) {
        t8_debugf ("Computing element for level %i, Max is %i\n", element_level, T8_DLINE_MAXLEVEL);
        // TODO: The maxlevel eventually should be element dependent. I know that an element dependent maxlevel function was developed in a different branch (by Sandro?)
        const t8_linearidx_t element_id_at_next_level
          = scheme->element_get_linear_id (tree_class, leaf_element, element_level + 1);

        const t8_locidx_t search_index
          = t8_forest_bin_search_lower (leafs, element_id_at_next_level, element_level + 1);
        // We expect the leaf element to be found at position ielement
        EXPECT_EQ (search_index, ielement)
          << "Found wrong position of level " << element_level + 1 << " leaf element with id "
          << element_id_at_next_level << ". Expected: " << ielement << " got: " << search_index;
      }

      // Construct an element that is definitely not in the array and
      // does not have an element of smaller id in the array. We expect -1 as return.
      // We take the first element of the forest and subtract 1 from its id.
      if (ielement == 0 && element_id > 0) {
        const t8_linearidx_t element_not_found_id = element_id - 1;
        const t8_locidx_t search_index = t8_forest_bin_search_lower (leafs, element_not_found_id, element_level);
        EXPECT_EQ (search_index, -1) << "Wrong return value for element that should not be in array. Expectec -1.";
      }
    }
  }
}

// bin_search_lower test for uniform forest
TEST_P (t8_bin_search_tester, bin_search_lower_uniform)
{
  t8_test_forest_bin_search_lower (forest);
}

// bin_search_lower test for adaptive forest
TEST_P (t8_bin_search_tester, bin_search_lower_adapt)
{
  t8_test_forest_bin_search_lower (forest_adapt);
}

// bin_search_upper test for uniform forest
TEST_P (t8_bin_search_tester, bin_search_upper_uniform)
{
  t8_test_forest_bin_search_upper (forest);
}

// bin_search_upper test for adaptive forest
TEST_P (t8_bin_search_tester, bin_search_upper_adapt)
{
  t8_test_forest_bin_search_upper (forest_adapt);
}

// TODO: Add these lambda to common headers since it is reused

/* Define a lambda to beatify gtest output for tuples <level, cmesh>.
 * This will set the correct level and cmesh name as part of the test case name. */
auto pretty_print_eclass_scheme_and_level
  = [] (const testing::TestParamInfo<std::tuple<std::tuple<int, t8_eclass_t>, int>> &info) {
      std::string scheme = t8_scheme_to_string[std::get<0> (std::get<0> (info.param))];
      std::string eclass = t8_eclass_to_string[std::get<1> (std::get<0> (info.param))];
      std::string level = std::string ("_level_") + std::to_string (std::get<1> (info.param));
      return scheme + "_" + eclass + level;
    };

INSTANTIATE_TEST_SUITE_P (t8_gtest_bin_search, t8_bin_search_tester,
                          testing::Combine (AllSchemes, testing::Range (0, T8_IS_LEAF_MAX_LVL)),
                          pretty_print_eclass_scheme_and_level);
