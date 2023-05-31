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
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include "t8_cmesh/t8_cmesh_testcases.h"
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <bitset>

#define MAX_NUM_ELEMETS 32

/* In this test, a uniform forest with x many elements is given. 
 * There are 2^x many ways to remove these x elements. After removing the elements,
 * the forest is coarsened and checked for overlapping elements.
 * 
 * Example:
 * num of elements = 4
 * instances - bianry with num elements many digits
 *      0    -    0 0 0 0
 *      1    -    0 0 0 1
 *      2    -    0 0 1 0
 *      3    -    0 0 1 1
 *      4    -    0 1 0 0
 *      5    -    0 1 0 1
 *      6    -    0 1 1 0
 *      7    -    0 1 1 1
 *      8    -    1 0 0 0
 *      9    -    1 0 0 1
 *      10   -    1 0 1 0
 *      11   -    1 0 1 1
 *      12   -    1 1 0 0
 *      13   -    1 1 0 1
 *      14   -    1 1 1 0
 *      15   -    1 1 1 1
 * We remove every element with id i if the i`th bit in the current instances is 0.
 */

/* *INDENT-OFF* */
class forest_permute:public testing::TestWithParam <t8_eclass_t> {
protected:
  void SetUp () override {
    eclass = GetParam();
#if T8_ENABLE_LESS_TESTS
    level = 1;
#else
    level = eclass < 4 ? 2 : 1;
#endif
    forest =
        t8_forest_new_uniform (t8_cmesh_new_from_class
                               (eclass, sc_MPI_COMM_WORLD), 
                               t8_scheme_new_default_cxx (),
                               level, 0, sc_MPI_COMM_WORLD);
    
    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &MPI_size);
    if (MPI_size > 1) {
      GTEST_SKIP ();
    } 
  }
  void TearDown () override {
    t8_forest_unref (&forest);
  }
  int                 MPI_size;
  int                 level;
  t8_eclass_t         eclass;
  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
};
/* *INDENT-ON* */

/** Remove every element with local_id i if the i`th bit in the current permutation \a removes is 0. */
static int
t8_adapt_remove (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts,
                 const int is_family,
                 const int num_elements,
                 t8_element_t *elements[])
{
  std::bitset<MAX_NUM_ELEMETS> *removes = (std::bitset<MAX_NUM_ELEMETS>*) t8_forest_get_user_data (forest);
  if (removes[lelement_id] == 0) {
    return -2;
  }
  return 0;
}

static int
t8_adapt_coarse (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t which_tree,
                 t8_locidx_t lelement_id,
                 t8_eclass_scheme_c *ts,
                 const int is_family,
                 const int num_elements,
                 t8_element_t *elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from,
                 t8_forest_adapt_t adapt_fn,
                 int do_partition, 
                 void *user_data)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, 0);
  
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  if (do_partition) {
    t8_forest_set_partition (forest_new, NULL, 0);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (forest_permute, test_permute_hole)
{
  /* number of instances */
  const t8_locidx_t   num_elments = t8_forest_get_tree_num_elements (forest, 0);
  T8_ASSERT (num_elments < MAX_NUM_ELEMETS);
  const uint32_t      num_permutation = 1 << num_elments;

  for (uint32_t permutation = 1; permutation < num_permutation; permutation++)
  {
    std::bitset<MAX_NUM_ELEMETS> removes (permutation);
    t8_forest_ref (forest);
    t8_forest_t         forest_adapt = t8_adapt_forest (forest, t8_adapt_remove, 0, &removes);
    for (int level = 0; level < 2; level++) {
      forest_adapt = t8_adapt_forest (forest_adapt, t8_adapt_coarse, 0, NULL);
    }
    ASSERT_TRUE (t8_forest_no_overlap (forest));
    t8_forest_unref (&forest_adapt);
  }
}

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_permute_hole, forest_permute, testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT));
/* *INDENT-ON* */