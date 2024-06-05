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
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>


/* In this test we check the t8_forest_element_is_leaf function.
 * Iterating over all cmesh test cases, we creat a uniform and an adaptive forest.
 * For each forest, we check that for each leaf element t8_forest_element_is_leaf returns true
 * and that it returns false for the parent and the first child.
 */


/* Maximum uniform level for forest. */
#define T8_IS_LEAF_MAX_LVL 4

/* Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_first_child (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                       t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  T8_ASSERT (!is_family || (is_family && num_elements == ts->t8_element_num_children (elements[0])));

  int level = ts->t8_element_level (elements[0]);

  /* we set a maximum refinement level as forest user data */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  int child_id = ts->t8_element_child_id (elements[0]);
  if (child_id == 1) {
    return 1;
  }
  return 0;
}


class element_is_leaf: public testing::TestWithParam<std::tuple<int,cmesh_example_base *>> {
 protected:

  void
  SetUp () override
  {
    /* Construct a cmesh */
    const int level = std::get<0>(GetParam ());
    cmesh = std::get<1>(GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes, we skip this case */
      GTEST_SKIP ();
    }
    /* Build the default scheme (TODO: Test this with all schemes) */
    scheme = t8_scheme_new_default_cxx ();
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
    if (t8_cmesh_is_empty (cmesh)) {
      t8_cmesh_destroy (&cmesh);
    }
    else {
      t8_forest_unref (&forest);
      t8_forest_unref (&forest_adapt);
    }
  }

  t8_cmesh_t cmesh;
  t8_forest_t forest, forest_adapt;
  t8_scheme_cxx_t *scheme;
};


void
t8_test_element_is_leaf_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  for (t8_locidx_t itree = 0;itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_eclass_scheme_c * scheme = t8_forest_get_eclass_scheme (forest, tree_class);
    /* Allocate memory to build a non-leaf element. */
    t8_element_t *not_leaf;
    scheme->t8_element_new (1, &not_leaf);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf returns false. */
    for (t8_locidx_t ielement = 0;ielement < num_elements_in_tree;++ielement) {
      const t8_element_t * leaf_element = t8_forest_get_element_in_tree (forest, itree, ielement);
      EXPECT_TRUE (t8_forest_element_is_leaf (forest, leaf_element, itree));
      /* Compute parent and first child of element and check that they are not in the tree */
      const int element_level = scheme->t8_element_level (leaf_element);
      if (element_level > 0) {
        scheme->t8_element_parent (leaf_element, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
      }
      if (element_level < scheme->t8_element_maxlevel ()) {
        scheme->t8_element_child (leaf_element, 0, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
      }
    }
    scheme->t8_element_destroy (1, &not_leaf);
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

/* TODO: Figure out how to do pretty printing for the cmesh examples when they are in a combined parameter. */
INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_leaf, element_is_leaf, testing::Combine (testing::Range (0, T8_IS_LEAF_MAX_LVL), AllCmeshsParam));
