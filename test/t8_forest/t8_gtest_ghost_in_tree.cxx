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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>
#include <t8_data/t8_element_array_iterator.hxx>

/* This test program tests the t8_ghost_get_ghost_in_tree_from_linear_id
 * and t8_ghost_get_ghost_id_in_tree function of the Ghost interface.
 * We adapt a forest and create its ghost layer. Afterwards, we
 * parse through all ghost elements and perform checks.
  */

static int
t8_test_adapt_second_element (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                              t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                              const int num_elements, t8_element_t *elements[])
{
  /* refine every second element up to the provided maximum level */
  int level = ts->t8_element_level (elements[0]);
  t8_linearidx_t eid = ts->t8_element_get_linear_id (elements[0], level);
  int maxlevel = *(int *) t8_forest_get_user_data (forest);

  if (eid % 2 && level < maxlevel) {
    return 1;
  }
  return 0;
}

class forest_ghost_in_tree: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {

    /* Construct a cmesh */
    cmesh = GetParam ()->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* empty cmeshes do not have ghosts and are skipped */
      GTEST_SKIP ();
    }

    scheme = t8_scheme_new_default_cxx ();
    /* Create a uniformly refined forest */
    t8_forest_t forest_uniform = t8_forest_new_uniform (cmesh, scheme, level, 1, sc_MPI_COMM_WORLD);
    /* Adapt the forest */
    int maxlevel = level + 2;
    forest_adapt = t8_forest_new_adapt (forest_uniform, t8_test_adapt_second_element, 1, 1, &maxlevel);
    // We need to update the cmesh since it changes with adapt.
    cmesh = t8_forest_get_cmesh (forest_adapt);
  }
  void
  TearDown () override
  {
    if (t8_cmesh_is_empty (cmesh)) {
      t8_cmesh_destroy (&cmesh);
    }
    else {
      t8_forest_unref (&forest_adapt);
    }
  }
  t8_cmesh_t cmesh;
  t8_scheme_cxx_t *scheme;
  t8_forest_t forest_adapt;
  int level =  // initial uniform refinement level
#if T8_ENABLE_LESS_TESTS
    3;  // less-tests, shorter runtime
#else
    5;  // more tests, longer runtime
#endif
};

// For reference:
// t8_ghost_get_ghost_id_in_tree
// Retrieves the local index of a ghost element in its specific ghost tree.
//
// t8_ghost_get_ghost_in_tree_from_linear_id
// Retrieves a ghost element from its ghost tree given the element's linear id.

TEST_P (forest_ghost_in_tree, test_get_ghost_id_in_tree)
{
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest_adapt);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest_adapt);

  for (t8_locidx_t ighost_tree = 0; ighost_tree < num_ghost_trees; ++ighost_tree) {
    // iterate over all Ghost trees
    // Get the ghost elements of that tree
    const t8_element_array_t *ghosts = t8_forest_ghost_get_tree_elements (forest_adapt, ighost_tree);
    t8_locidx_t ghost_index_in_tree = 0;
    //for (const t8_element_t &ighost : ghosts) {  // loop over all ghost elements

    for (; ghost_index_in_tree < t8_element_array_get_count (ghosts); ++ghost_index_in_tree) {
      const t8_element_t *ighost = t8_element_array_index_locidx (ghosts, ghost_index_in_tree);
      // Compute the index in the tree of the current ghost
      const t8_locidx_t check_ghost_index_in_tree = t8_ghost_get_ghost_id_in_tree (forest_adapt, ighost_tree, ighost);
      ASSERT_GE (check_ghost_index_in_tree, 0) << "Ghost element was not found.";  // The ghost must have been found.
      // Check that the index is correct
      EXPECT_EQ (check_ghost_index_in_tree, ghost_index_in_tree) << "Wrong index returned for ghost element.";
      // ++ghost_index_in_tree;
    }
  }
}

TEST_P (forest_ghost_in_tree, test_get_ghost_in_tree_from_linear_id)
{
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest_adapt);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest_adapt);

  for (t8_locidx_t ighost_tree = 0; ighost_tree < num_ghost_trees; ++ighost_tree) {
    // iterate over all Ghost trees
    // Get the ghost elements of that tree
    const t8_element_array_t *ghosts = t8_forest_ghost_get_tree_elements (forest_adapt, ighost_tree);
    t8_locidx_t ghost_index_in_tree = 0;
    //for (const t8_element_t &ighost : ghosts) {  // loop over all ghost elements
    for (; ghost_index_in_tree < t8_element_array_get_count (ghosts); ++ghost_index_in_tree) {

      const t8_element_t *ighost = t8_element_array_index_locidx (ghosts, ghost_index_in_tree);
      // Now checking t8_ghost_get_ghost_in_tree_from_linear_id
      const t8_eclass_t eclass = t8_forest_get_tree_class (forest_adapt, num_local_trees + ighost_tree);
      const t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest_adapt, eclass);
      const int ghost_level = scheme->t8_element_level (ighost);
      const t8_linearidx_t ghost_id = scheme->t8_element_get_linear_id (ighost, ghost_level);
      t8_locidx_t check_ghost_index;
      const t8_element_t *check_ghost = t8_ghost_get_ghost_in_tree_from_linear_id (forest_adapt, ighost_tree, ghost_id,
                                                                                   ghost_level, &check_ghost_index);

      EXPECT_TRUE (scheme->t8_element_equal (ighost, check_ghost)) << "Returned wrong ghost element.";
      EXPECT_EQ (check_ghost_index, ghost_index_in_tree) << "Returned wrong ghost index.";

      // ++ghost_index_in_tree;
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_ghost_in_tree, forest_ghost_in_tree, AllCmeshsParam,
                          pretty_print_base_example);
