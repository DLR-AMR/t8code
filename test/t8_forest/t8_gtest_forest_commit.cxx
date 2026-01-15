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
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

/* In this test we adapt, balance and partition a uniform forest.
 * We do this in two ways:
 * 1st  All operations are performed in one single call to t8_forest_commit
 * 2nd  Each intermediate step is performed in a separate commit
 *
 * After these two forests are created, we check for equality.
 */

struct forest_commit: public testing::TestWithParam<std::tuple<int, cmesh_example_base *>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    /* Construct a cmesh */
    cmesh = std::get<1> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes*/
      scheme->unref ();
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }
  t8_cmesh_t cmesh;
  const t8_scheme *scheme;
};

/* Adapt a forest such that always the first child of a
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_balance (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] t8_locidx_t which_tree, t8_eclass_t tree_class,
                       [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                       [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                       t8_element_t *elements[])
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
#if T8_TEST_LEVEL_INT < 2
  t8_forest_set_balance (forest_ada_bal_par, NULL, 0);
#endif
  t8_forest_set_partition (forest_ada_bal_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_par);

  return forest_ada_bal_par;
}

/* adapt, balance and partition a given forest in 3 steps */
static t8_forest_t
t8_test_forest_commit_abp_3step (t8_forest_t forest, int maxlevel)
{
  t8_forest_t forest_adapt;
  t8_forest_t forest_partition;
#if T8_TEST_LEVEL_INT < 2
  t8_forest_t forest_balance;
  t8_forest_init (&forest_balance);
#endif

  t8_forest_init (&forest_adapt);
  t8_forest_init (&forest_partition);

  /* adapt the forest */
  t8_forest_set_user_data (forest_adapt, &maxlevel);
  t8_forest_set_adapt (forest_adapt, forest, t8_test_adapt_balance, 1);
  t8_forest_commit (forest_adapt);

#if T8_TEST_LEVEL_INT < 2
  /* balance the forest */
  t8_forest_set_balance (forest_balance, forest_adapt, 0);
  t8_forest_commit (forest_balance);
#endif

  /* partition the forest */
#if T8_TEST_LEVEL_INT < 2
  t8_forest_set_partition (forest_partition, forest_balance, 0);
#else
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
#endif
  t8_forest_commit (forest_partition);

  return forest_partition;
}

TEST_P (forest_commit, test_forest_commit)
{

  t8_forest_t forest;
  t8_forest_t forest_ada_bal_part;
  t8_forest_t forest_abp_3part;

  const int level_step = 2;

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
    scheme->ref ();
    t8_forest_unref (&forest_ada_bal_part);
    t8_forest_unref (&forest_abp_3part);
  }
  scheme->unref ();
  t8_debugf ("Done testing forest commit.");
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_commit, forest_commit,
                          testing::Combine (AllSchemeCollections, AllCmeshsParam), pretty_print_base_example_scheme);
