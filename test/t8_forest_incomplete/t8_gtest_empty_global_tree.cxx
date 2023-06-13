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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* *INDENT-OFF* */
class global_tree:public testing::TestWithParam <t8_eclass_t> {
protected:
  void SetUp () override {
    eclass = GetParam();
    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &MPI_size);
    sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &MPI_rank);
    
    forest = 
      t8_forest_new_uniform (t8_cmesh_new_bigmesh (eclass, 3, sc_MPI_COMM_WORLD), t8_scheme_new_default_cxx (), 0, 0, sc_MPI_COMM_WORLD);  
  }
  void TearDown () override {
    t8_forest_unref (&forest);
  }
  int                 MPI_size;
  int                 MPI_rank;
  t8_eclass_t         eclass;
  t8_forest_t         forest;
};
/* *INDENT-ON* */

/**  */
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
  const int *instance = (const int *) t8_forest_get_user_data (forest);
  const t8_gloidx_t global_tree_id = t8_forest_global_tree_id (forest_from, which_tree);
  switch (*instance)
  {
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
    return -2;
    break;
  }
  
  return 0;
}

t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 void *user_data)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, 0);
  if (user_data != NULL) {
    t8_forest_set_user_data (forest_new, user_data);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

TEST_P (global_tree, test_empty_global_tree)
{
  ASSERT_TRUE (!forest->incomplete_trees);
  for (int instance = 0; instance < 7; instance++) {
    t8_forest_ref (forest);
    t8_forest_t         forest_adapt = t8_adapt_forest (forest, t8_adapt_remove, &instance);
    ASSERT_TRUE (forest_adapt->incomplete_trees);
    ASSERT_TRUE (!forest->incomplete_trees);
    t8_forest_unref (&forest_adapt);
  }
}

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_empty_global_tree, global_tree, testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT));
/* *INDENT-ON* */
