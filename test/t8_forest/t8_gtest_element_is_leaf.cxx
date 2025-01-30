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
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>

/* In this test we check the t8_forest_element_is_leaf function.
 * Iterating over all cmesh test cases, we creat a uniform and an adaptive forest.
 * For each forest, we check that for each leaf element t8_forest_element_is_leaf returns true
 * and that it returns false for the parent and the first child.
 */

/* Maximum uniform level for forest. */

#ifdef T8_ENABLE_LESS_TESTS
#define T8_IS_LEAF_MAX_LVL 3
#else
#define T8_IS_LEAF_MAX_LVL 4
#endif
/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_first_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                           const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme *scheme,
                           const int is_family, const int num_elements, t8_element_t *elements[])
{
  T8_ASSERT (!is_family || (is_family && num_elements == scheme->element_get_num_children (tree_class, elements[0])));

  const int level = scheme->element_get_level (tree_class, elements[0]);

  /* we set a maximum refinement level as forest user data */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  const int child_id = scheme->element_get_child_id (tree_class, elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}

class element_is_leaf: public testing::TestWithParam<std::tuple<int, int, cmesh_example_base *>> {
 protected:
  void
  SetUp () override
  {
    /* Construct a cmesh */
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    const int level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = std::get<2> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes, we skip this case */
      scheme->unref ();
      t8_cmesh_unref (&cmesh);
      GTEST_SKIP ();
    }
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

void
t8_test_element_is_leaf_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Allocate memory to build a non-leaf element. */
    t8_element_t *not_leaf;
    scheme->element_new (tree_class, 1, &not_leaf);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_element_in_tree (forest, itree, ielement);
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

/* Define a lambda to beatify gtest output for tuples <level, cmesh>.
 * This will set the correct level and cmesh name as part of the test case name. */
auto pretty_print_level_and_cmesh_params
  = [] (const testing::TestParamInfo<std::tuple<int, int, cmesh_example_base *>> &info) {
      std::string name = std::string ("Level_") + std::to_string (std::get<1> (info.param));
      std::string cmesh_name;
      std::get<2> (info.param)->param_to_string (cmesh_name);
      name += std::string ("_") + cmesh_name;
      name += std::string ("scheme_") + std::to_string (std::get<0> (info.param));
      name += std::string ("_") + std::to_string (info.index);
      return name;
    };

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_leaf, element_is_leaf,
                          testing::Combine (AllSchemeCollections, testing::Range (0, T8_IS_LEAF_MAX_LVL),
                                            AllCmeshsParam),
                          pretty_print_level_and_cmesh_params);
