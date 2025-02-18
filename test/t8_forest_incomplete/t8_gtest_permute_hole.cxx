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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <bitset>
#include <test/t8_gtest_macros.hxx>

#define MAX_NUM_ELEMENTS 32 /* number of digits for binary representation */

/* In this test, a uniform forest with x many elements is given. 
 * There are 2^x many ways to remove these x elements. After removing the elements,
 * the forest is coarsened and checked for overlapping elements.
 * 
 * Example:
 * num of elements = 4
 * instances - binary representation
 *      0    -  ...0 0 0 0 0 0 0
 *      1    -  ...0 0 0 0 0 0 1
 *      2    -  ...0 0 0 0 0 1 0
 *      3    -  ...0 0 0 0 0 1 1
 *      4    -  ...0 0 0 0 1 0 0
 *      5    -  ...0 0 0 0 1 0 1
 *      6    -  ...0 0 0 0 1 1 0
 *      7    -  ...0 0 0 0 1 1 1
 *      8    -  ...0 0 0 1 0 0 0
 *      9    -  ...0 0 0 1 0 0 1
 *      10   -  ...0 0 0 1 0 1 0
 *      11   -  ...0 0 0 1 0 1 1
 *      12   -  ...0 0 0 1 1 0 0
 *      13   -  ...0 0 0 1 1 0 1
 *      14   -  ...0 0 0 1 1 1 0
 *      15   -  ...0 0 0 1 1 1 1
 * We remove every element with id i if the i`th bit (right to left) 
 * in the current instances is 0.
 * 
 * Note, this test runs only on one rank.
 */

class forest_permute: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
#if T8_ENABLE_LESS_TESTS
    level = 1;
#else
    level = eclass < 4 ? 2 : 1;
#endif
    forest = t8_forest_new_uniform (t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD), t8_scheme_new_default (),
                                    level, 0, sc_MPI_COMM_WORLD);

    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &MPI_size);
    if (MPI_size > 1) {
      GTEST_SKIP ();
    }
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  int MPI_size;
  int level;
  t8_eclass_t eclass;
  t8_forest_t forest;
  t8_cmesh_t cmesh;
};

/** This structure contains a bitset with all 
 *  elements to be removed.
 */
struct t8_elements
{
  std::bitset<MAX_NUM_ELEMENTS> remove;
};

/** Remove every element with local_id i if the i`th bit in 
 * the current permutation \a remove is 0. */
static int
t8_adapt_remove (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
                 t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme, [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                 [[maybe_unused]] t8_element_t *elements[])
{
  struct t8_elements *data = (struct t8_elements *) t8_forest_get_user_data (forest);
  if (data->remove[lelement_id] == 0) {
    return -2;
  }
  return 0;
}

/** Coarse every (incomplete) family */
static int
t8_adapt_coarse ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, [[maybe_unused]] t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
    [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme, const int is_family, [[maybe_unused]] const int num_elements,
    [[maybe_unused]] t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, void *user_data)
{
  t8_forest_t forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, 0);
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (forest_permute, test_permute_hole)
{
  /* number of instances/permutations */
  const t8_locidx_t num_elements = t8_forest_get_tree_num_elements (forest, 0);
  T8_ASSERT (num_elements < MAX_NUM_ELEMENTS);
  const uint32_t num_permutation = 1 << num_elements;

  for (uint32_t permutation = 1; permutation < num_permutation; permutation++) {
    std::bitset<MAX_NUM_ELEMENTS> remove (permutation);

    struct t8_elements data
    {
      remove
    };

    t8_forest_ref (forest);
    /* Remove elements for every 0 bit in \a removes */
    t8_forest_t forest_adapt = t8_adapt_forest (forest, t8_adapt_remove, &data);

    /* check if the correct number of elements got removed */
    t8_locidx_t element_count = 0;
    for (t8_locidx_t ridx = 0; ridx < MAX_NUM_ELEMENTS; ridx++) {
      if (remove[ridx] == 1) {
        element_count++;
      }
    }
    ASSERT_TRUE (element_count == t8_forest_get_tree_num_elements (forest_adapt, 0));

    /* check if coarsening results in overlapping elements */
    for (int l = 0; l < level + 1; l++) {
      forest_adapt = t8_adapt_forest (forest_adapt, t8_adapt_coarse, NULL);
      ASSERT_TRUE (t8_forest_no_overlap (forest_adapt));
    }
    ASSERT_TRUE (1 == t8_forest_get_tree_num_elements (forest_adapt, 0));

    t8_forest_unref (&forest_adapt);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_permute_hole, forest_permute, AllEclasses, print_eclass);
