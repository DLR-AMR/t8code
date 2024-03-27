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
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <sc_functions.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

class linear_id: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    scheme = t8_scheme_new_default_cxx ();
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test);
    ts->t8_element_root (element);
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &test);
    t8_scheme_cxx_unref (&scheme);
  }
  t8_element_t *element;
  t8_element_t *child;
  t8_element_t *test;
  t8_scheme_cxx *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

static int
t8_test_init_linear_id_refine_everything (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                          t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                                          const int num_elements, t8_element_t *elements[])
{
  return 1;
}

/* Iterate over the leaves of a uniformly refined forest and check the id*/
TEST_P (linear_id, uniform_forest)
{
  t8_forest_t forest, forest_adapt;
  t8_cmesh_t cmesh;
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 5;
#else
  const int maxlvl = 6;
#endif
  /* Construct a forest with a single element of the current class*/
  cmesh = t8_cmesh_new_from_class (ts->eclass, comm);
  t8_cmesh_ref (cmesh);
  forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  t8_scheme_cxx_ref (scheme);
  for (int level = 0; level < maxlvl; level++) {
    /*Get the number of local trees*/
    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
    /*Iterate over trees */
    for (t8_locidx_t tree_id = 0; tree_id < num_local_trees; tree_id++) {
      /*Get the number of elements in the tree*/
      const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, tree_id);
      /*Manually compute the id of the first element*/
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, tree_id);
      t8_eclass_scheme_c *tc_scheme = t8_forest_get_eclass_scheme (forest, tree_class);
      const t8_locidx_t shift = tc_scheme->t8_element_count_leaves_from_root (level) - num_elements_in_tree;
      /*Iterate over elements */
      for (t8_locidx_t id_iter = 0; id_iter < num_elements_in_tree; id_iter++) {
        /*Get the current element*/
        const t8_element_t *element = t8_forest_get_element_in_tree (forest, tree_id, id_iter);
        /*Get the ID of the element at current level */
        const t8_locidx_t id = ts->t8_element_get_linear_id (element, level);
        /* Check the computed id*/
        EXPECT_EQ (id, id_iter + shift);
      }
    }
    /* Construct the uniformly refined forest of the next level */
    t8_forest_init (&forest_adapt);
    t8_forest_set_level (forest_adapt, level + 1);
    t8_forest_set_adapt (forest_adapt, forest, t8_test_init_linear_id_refine_everything, 0);
    t8_forest_commit (forest_adapt);
    forest = forest_adapt;
  }
  t8_cmesh_unref (&cmesh);
  t8_forest_unref (&forest_adapt);
}

/* Test, if the linear_id of descendants of an element is the same as the id of element 
 * (on the level defined by the element) */
TEST_P (linear_id, id_at_other_level)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int max_lvl = 3; /* Maximal level to compute elements on */
  const int add_lvl = 3; /* maxlvl + add_lvl is the level of the descendants*/
#else
  const int max_lvl = 4;
  const int add_lvl = 3;
#endif
  for (int level = 0; level < max_lvl; level++) {
    /* Compute the number of elements at the current level */
    const t8_linearidx_t num_desc = ts->t8_element_count_leaves_from_root (level);
    for (t8_linearidx_t id = 0; id < num_desc; id++) {
      /* Set the child at the current level */
      ts->t8_element_set_linear_id (child, level, id);
      /* Compute the id of child at a higher level. */
      const t8_linearidx_t id_at_lvl = ts->t8_element_get_linear_id (child, level + add_lvl);
      /* Compute how many leaves/descendants child has at level level+add_lvl */
      const t8_linearidx_t child_desc = ts->t8_element_count_leaves (child, level + add_lvl);
      /* Iterate over all descendants */
      for (t8_linearidx_t leaf_id = 0; leaf_id < child_desc; leaf_id++) {
        /* Set the descendant (test) at level of the descendants and shift the 
         * leaf_id into the region of the descendants of child*/
        ts->t8_element_set_linear_id (test, level + add_lvl, id_at_lvl + leaf_id);
        /* Compute the id of the descendant (test) at the current level */
        const t8_linearidx_t test_id = ts->t8_element_get_linear_id (test, level);
        /* test_id and id should be equal. */
        EXPECT_EQ (id, test_id);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_test_init_linear_id, linear_id, AllEclasses, print_eclass);
