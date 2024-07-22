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
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>

/* In this test, we adapt, balance and partition a uniform forest.
 * We do this in two ways:
 * 1st  All operations are performed in one single call to t8_forest_commit
 * 2nd  Each intermediate step is performed in a separate commit
 *
 * After these two forests are created, we check for equality.
 */

class forest_commit: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    /* Construct a cmesh */
    cmesh = GetParam ()->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes*/
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }
  t8_cmesh_t cmesh;
};

/* Adapt a forest such that always the first child of a
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_balance (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
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

/* adapt, balance and partition a given forest in one step */
static t8_forest_t
t8_test_forest_commit_abp (t8_forest_t forest, int maxlevel)
{
  t8_forest_t forest_ada_bal_par;

  /* Adapt, balance and partition the uniform forest */
  t8_forest_init (&forest_ada_bal_par);
  /* Set user data for adapt */
  t8_forest_set_user_data (forest_ada_bal_par, &maxlevel);
  t8_forest_set_adapt (forest_ada_bal_par, forest, t8_test_adapt_balance, 1);
  t8_forest_set_balance (forest_ada_bal_par, NULL, 0);
  t8_forest_set_partition (forest_ada_bal_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_par);

  return forest_ada_bal_par;
}

/* adapt, balance and partition a given forest in 3 steps */
static t8_forest_t
t8_test_forest_commit_abp_3step (t8_forest_t forest, int maxlevel)
{
  t8_forest_t forest_adapt;
  t8_forest_t forest_balance;
  t8_forest_t forest_partition;

  t8_forest_init (&forest_adapt);
  t8_forest_init (&forest_balance);
  t8_forest_init (&forest_partition);

  /* adapt the forest */
  t8_forest_set_user_data (forest_adapt, &maxlevel);
  t8_forest_set_adapt (forest_adapt, forest, t8_test_adapt_balance, 1);
  t8_forest_commit (forest_adapt);

  /* balance the forest */
  t8_forest_set_balance (forest_balance, forest_adapt, 0);
  t8_forest_commit (forest_balance);

  /* partrition the forest */
  t8_forest_set_partition (forest_partition, forest_balance, 0);
  t8_forest_commit (forest_partition);

  return forest_partition;
}

TEST_P (forest_commit, test_forest_commit)
{

  t8_forest_t forest;
  t8_forest_t forest_ada_bal_part;
  t8_forest_t forest_abp_3part;

  const int level_step = 2;

  t8_scheme_cxx_t *scheme = t8_scheme_new_standalone_cxx ();

  /* Compute the first level, such that no process is empty */
  int min_level = t8_forest_min_nonempty_level (cmesh, scheme);
  /* Use one level with empty processes */
  min_level = SC_MAX (min_level - 1, 0);
  for (int level = min_level; level < min_level + level_step; level++) {
    t8_debugf ("Testing forest commit level %i\n", level);
    int maxlevel = level + level_step;
    /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    /* Create a uniformly refined forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 1, sc_MPI_COMM_WORLD);
    /* We need to use forest twice, so we ref it */
    t8_forest_ref (forest);
    /* Adapt, balance and partition the forest */
    forest_ada_bal_part = t8_test_forest_commit_abp (forest, maxlevel);
    /* Adapt, balance and partition the forest using three separate steps */
    forest_abp_3part = t8_test_forest_commit_abp_3step (forest, maxlevel);

    ASSERT_TRUE (t8_forest_is_equal (forest_abp_3part, forest_ada_bal_part)) << "The forests are not equal";
    t8_scheme_cxx_ref (scheme);
    t8_forest_unref (&forest_ada_bal_part);
    t8_forest_unref (&forest_abp_3part);
  }
  t8_scheme_cxx_unref (&scheme);
  t8_debugf ("Done testing forest commit.");
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_commit, forest_commit, AllCmeshsParam, pretty_print_base_example);
