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
#include <t8.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

/** In this test, we are given a forest with 3 global trees. 
 * We adapt the forest so that all 6 compositions of empty 
 * global trees are the result of it. 
 * Therefore, \a testcase runs from 0 to 5.
 * We do this twice. Once we partition the forest in the same call. 
 * The second time, we do the adapting and partitioning separately.
 * The two resulting forests must be equal.
 * */

/* Remove `DISABLED_` from the name of the Test(suite) or use `--gtest_also_run_disabled_tests` when you start working on the issue. */
class DISABLED_global_tree: public testing::TestWithParam<std::tuple<t8_eclass, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    testcase = std::get<1> (GetParam ());
    forest = t8_forest_new_uniform (t8_cmesh_new_bigmesh (eclass, 3, sc_MPI_COMM_WORLD), t8_scheme_new_default (), 0, 0,
                                    sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  int testcase;
  t8_eclass_t eclass;
  t8_forest_t forest;
};

/** Removes all elements of local trees if they belong to the corresponding
 *  global trees which are defined by the current testcase of test. */
static int
t8_adapt_remove (t8_forest_t forest_from, t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
                 [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme,
                 [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                 [[maybe_unused]] t8_element_t *elements[], [[maybe_unused]] void *user_data)
{
  const int *testcase = (const int *) user_data;
  const t8_gloidx_t global_tree_id = t8_forest_global_tree_id (forest_from, which_tree);
  switch (*testcase) {
  case 0:
    if (global_tree_id == 0) {
      return -2;
    }
    break;
  case 1:
    if (global_tree_id == 1) {
      return -2;
    }
    break;
  case 2:
    if (global_tree_id == 2) {
      return -2;
    }
    break;
  case 3:
    if (global_tree_id != 0) {
      return -2;
    }
    break;
  case 4:
    if (global_tree_id != 1) {
      return -2;
    }
    break;
  case 5:
    if (global_tree_id != 2) {
      return -2;
    }
    break;
  default:
    break;
  }

  return 0;
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int do_adapt, int do_partition, void *user_data)
{
  t8_forest_t forest_new;

  t8_forest_init (&forest_new);
  if (do_adapt) {
    t8_forest_set_adapt (forest_new, forest_from, adapt_fn, 0);
    if (do_partition) {
      t8_forest_set_partition (forest_new, NULL, 0);
    }
  }
  else if (do_partition) {
    t8_forest_set_partition (forest_new, forest_from, 0);
  }
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (DISABLED_global_tree, test_empty_global_tree)
{
  ASSERT_TRUE (!forest->incomplete_trees);

  t8_forest_ref (forest);
  /* Do adapt and partition in one step */
  t8_forest_t forest_adapt_a = t8_adapt_forest (forest, t8_adapt_remove, 1, 1, &testcase);
  ASSERT_TRUE (forest_adapt_a->incomplete_trees);
  ASSERT_TRUE (!forest->incomplete_trees);

  t8_forest_ref (forest);
  /* Do adapt and partition in separate steps */
  t8_forest_t forest_adapt_b = t8_adapt_forest (forest, t8_adapt_remove, 1, 0, &testcase);
  ASSERT_TRUE (forest_adapt_b->incomplete_trees);
  ASSERT_TRUE (!forest->incomplete_trees);
  forest_adapt_b = t8_adapt_forest (forest_adapt_b, NULL, 0, 1, NULL);

  /* The number of trees and elements between forest_adapt_a and forest_adapt_b have to match. */
  /* Global */
  ASSERT_EQ (t8_forest_get_num_global_trees (forest), t8_forest_get_num_global_trees (forest_adapt_a));
  ASSERT_EQ (t8_forest_get_num_global_trees (forest_adapt_b), t8_forest_get_num_global_trees (forest_adapt_a));
  ASSERT_EQ (t8_forest_get_global_num_leaf_elements (forest_adapt_b),
             t8_forest_get_global_num_leaf_elements (forest_adapt_a));

  /* Compare forest->global_num_trees with the sum of all local trees
   * on all processes. Those numbers must be equal, since every tree
   * contains only the root element and no tree is shared. */
  t8_gloidx_t local_num_trees = (t8_gloidx_t) t8_forest_get_num_local_trees (forest_adapt_a);
  t8_gloidx_t global_num_trees;
  int mpiret
    = sc_MPI_Allreduce (&local_num_trees, &global_num_trees, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM, sc_MPI_COMM_WORLD);
  SC_CHECK_MPI (mpiret);
  ASSERT_EQ (global_num_trees, t8_forest_get_num_global_trees (forest));

  /* Local */
  ASSERT_TRUE (t8_forest_is_equal (forest_adapt_b, forest_adapt_a));

  t8_forest_unref (&forest_adapt_a);
  t8_forest_unref (&forest_adapt_b);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_empty_global_tree, DISABLED_global_tree,
                          testing::Combine (AllEclasses, testing::Range (0, 6)));
