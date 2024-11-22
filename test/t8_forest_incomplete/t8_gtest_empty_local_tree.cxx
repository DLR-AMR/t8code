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
#include <bitset>

#define MAX_NUM_RANKS 8

/* In this test, a partitioned forest with one global tree and at 
 * least so many elements, such that each process has at least one 
 * local element is given. Let x be the number of mpi ranks.
 * There are 2^x many ways to empty these x local trees.
 * 
 * Example:
 * x = 3
 * instances - binary representation
 *      0    -  0 0
 *      1    -  0 1
 *      2    -  1 0
 *      3    -  1 1
 * We remove all elements from rank with id i if the i`th bit 
 * in the current instances is 0.
 * 
 * Note, this test runs only on two to maxmal 8 ranks.
 * 
 * We adapt the given forest twice. 
 * The first time, we partition the forest in the same call. 
 * The second time, we do the adapting and partitioning separately.
 * The two resulting forests must be equal.
 */

/** This test covers the functionality described in Issue: https://github.com/DLR-AMR/t8code/issues/1137
 * Remove `DISABLED_` from the name of the Test(suite) or use `--gtest_also_run_disabled_tests` when you start working on the issue. */
class DISABLED_local_tree: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &MPI_size);

    forest = t8_forest_new_uniform (t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD), t8_scheme_new_default (),
                                    MPI_size, 0, sc_MPI_COMM_WORLD);
    /* TODO: The level does not need to be as big as MPI_SIZE, only as big so that each process has at least one element */

    if (MPI_size == 1 || MPI_size > MAX_NUM_RANKS) {
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  int MPI_size;
  t8_eclass_t eclass;
  t8_forest_t forest;
};

/** This structure contains a bitset with all 
 * local trees on all processes to be removed.
 */
struct t8_trees_to_remove
{
  std::bitset<MAX_NUM_RANKS> remove;
};

/** Remove every element of rank i if the i`th bit in 
 * the current instance \a remove is 0. */
static int
t8_adapt_remove (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
                 t8_locidx_t lelement_id, const t8_scheme *ts, const int is_family, const int num_elements,
                 t8_element_t *elements[])
{
  struct t8_trees_to_remove *trees_to_remove = (struct t8_trees_to_remove *) t8_forest_get_user_data (forest);
  if (trees_to_remove->remove[forest_from->mpirank] == 0) {
    return -2;
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

TEST_P (DISABLED_local_tree, test_empty_local_tree)
{
  /* Number of instances/testcases */
  const uint32_t num_instances = 1 << MPI_size;
  for (uint32_t instances = 1; instances < num_instances - 1; instances++) {
    /* Remove all local elements for every process \a remove[rank] == 0 */
    std::bitset<MAX_NUM_RANKS> remove (instances);

    struct t8_trees_to_remove data
    {
      remove
    };

    t8_forest_ref (forest);
    /* Do adapt and partition in one step */
    t8_forest_t forest_adapt_a = t8_adapt_forest (forest, t8_adapt_remove, 1, 1, &data);
    ASSERT_TRUE (forest_adapt_a->incomplete_trees);
    ASSERT_TRUE (!forest->incomplete_trees);

    t8_forest_ref (forest);
    /* Do adapt and partition in separate steps */
    t8_forest_t forest_adapt_b = t8_adapt_forest (forest, t8_adapt_remove, 1, 0, &data);
    ASSERT_TRUE (forest_adapt_b->incomplete_trees);
    ASSERT_TRUE (!forest->incomplete_trees);
    forest_adapt_b = t8_adapt_forest (forest_adapt_b, NULL, 0, 1, NULL);

    /* The trees have to be equal */
    ASSERT_TRUE (t8_forest_is_equal (forest_adapt_b, forest_adapt_a));

    t8_forest_unref (&forest_adapt_a);
    t8_forest_unref (&forest_adapt_b);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_empty_local_tree, DISABLED_local_tree,
                          testing::Range (T8_ECLASS_LINE, T8_ECLASS_COUNT), print_eclass);
