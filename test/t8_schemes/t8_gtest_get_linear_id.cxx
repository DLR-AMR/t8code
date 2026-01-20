/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <t8_gtest_schemes.hxx>
#include <sc_functions.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_gtest_macros.hxx>

class get_linear_id: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    scheme->element_new (eclass, 1, &element);
    scheme->element_new (eclass, 1, &child);
    scheme->element_new (eclass, 1, &test);
    scheme->set_to_root (eclass, element);
  }

  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &element);
    scheme->element_destroy (eclass, 1, &child);
    scheme->element_destroy (eclass, 1, &test);
    scheme->unref ();
  }
  t8_element_t *element;
  t8_element_t *child;
  t8_element_t *test;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

static int
t8_test_init_linear_id_refine_everything ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                                          [[maybe_unused]] t8_locidx_t which_tree,
                                          [[maybe_unused]] const t8_eclass_t tree_class,
                                          [[maybe_unused]] t8_locidx_t lelement_id,
                                          [[maybe_unused]] const t8_scheme *scheme,
                                          [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                                          [[maybe_unused]] t8_element_t *elements[])
{
  return 1;
}

/* Iterate over the leaves of a uniformly refined forest and check the id*/
TEST_P (get_linear_id, uniform_forest)
{
  t8_forest_t forest, forest_adapt;
  t8_cmesh_t cmesh;
#if T8_TEST_LEVEL_INT >= 2
  const int maxlvl = 4;
#elif T8_TEST_LEVEL_INT >= 1
  const int maxlvl = 5;
#else
  const int maxlvl = 6;
#endif
  /* Construct a forest with a single element of the current class*/
  cmesh = t8_cmesh_new_from_class (eclass, comm);
  t8_cmesh_ref (cmesh);
  forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  const t8_scheme *tc_scheme = t8_forest_get_scheme (forest);
  scheme->ref ();
  for (int level = 0; level < maxlvl; level++) {
    /*Get the number of local trees*/
    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
    /*Iterate over trees */
    for (t8_locidx_t tree_id = 0; tree_id < num_local_trees; tree_id++) {
      /*Get the number of elements in the tree*/
      const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, tree_id);
      /*Manually compute the id of the first element*/
      const t8_linearidx_t shift = tc_scheme->count_leaves_from_root (eclass, level) - num_elements_in_tree;
      /*Iterate over elements */
      for (t8_locidx_t id_iter = 0; id_iter < num_elements_in_tree; id_iter++) {
        /*Get the current element*/
        const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, tree_id, id_iter);
        /*Get the ID of the element at current level */
        const t8_element_level elem_level = tc_scheme->element_get_level (eclass, element);
        const t8_linearidx_t id = tc_scheme->element_get_linear_id (eclass, element, elem_level);
        /* Check the computed id*/
        EXPECT_EQ (id, static_cast<t8_linearidx_t> (id_iter) + shift);
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
TEST_P (get_linear_id, id_at_other_level)
{
#if T8_TEST_LEVEL_INT >= 1
  const int max_lvl = 3; /* Maximal level to compute elements on */
  const int add_lvl = 3; /* maxlvl + add_lvl is the level of the descendants*/
#else
  const int max_lvl = 4;
  const int add_lvl = 3;
#endif
  for (int level = 0; level < max_lvl; level++) {
    /* Compute the number of elements at the current level */
    const t8_linearidx_t num_desc = scheme->count_leaves_from_root (eclass, level);
    for (t8_linearidx_t id = 0; id < num_desc; id++) {
      /* Set the child at the current level */
      scheme->element_set_linear_id (eclass, child, level, id);
      /* Compute the id of child at a higher level. */
      const t8_linearidx_t id_at_lvl = scheme->element_get_linear_id (eclass, child, level + add_lvl);
      /* Compute how many leaves/descendants child has at level level+add_lvl */
      const t8_linearidx_t child_desc = scheme->element_count_leaves (eclass, child, level + add_lvl);
      /* Iterate over all descendants */
      for (t8_linearidx_t leaf_id = 0; leaf_id < child_desc; leaf_id++) {
        /* Set the descendant (test) at level of the descendants and shift the
         * leaf_id into the region of the descendants of child*/
        scheme->element_set_linear_id (eclass, test, level + add_lvl, id_at_lvl + leaf_id);
        /* Compute the id of the descendant (test) at the current level */
        const t8_linearidx_t test_id = scheme->element_get_linear_id (eclass, test, level);
        /* test_id and id should be equal. */
        EXPECT_EQ (id, test_id);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_test_get_linear_id, get_linear_id, AllSchemes, print_all_schemes);
