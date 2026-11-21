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
#include <t8_gtest_schemes.hxx>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include <t8_gtest_adapt_callbacks.hxx>
#include <t8_gtest_macros.hxx>

/** \file t8_gtest_element_is_leaf.cxx
 *
 * For each forest, we check that for each leaf (ghost) element t8_forest_element_is_leaf
 * (t8_forest_element_is_ghost) returns true and that it returns false for the parent and the first child.
 *
 * Fixture element_is_leaf: Since t8_forest_element_is_leaf works on a single tree, we do not need
 * to test all example cmeshes. For each scheme and element class, we create a single tree cmesh
 * and build a uniform and an adaptive forest without ghost layer on it.
 *
 * Fixture element_is_ghost: Ghost elements depend on multiple trees and processes. Iterating over
 * all cmesh test cases, we create a uniform and an adaptive forest, both with a ghost layer.
 */

/* Maximum uniform level for forest. */

#if T8_TEST_LEVEL_INT >= 1
#define T8_IS_LEAF_MAX_LVL 2
#else
#define T8_IS_LEAF_MAX_LVL 3
#endif

struct element_is_leaf: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    const t8_eclass_t eclass = std::get<1> (std::get<0> (GetParam ()));
    const int level = std::get<1> (GetParam ());
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);

    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_class (cmesh, eclass, sc_MPI_COMM_WORLD);

    const int do_ghost = 0;
    forest = t8_forest_new_uniform (cmesh, scheme, level, do_ghost, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    int maxlevel = level + 2;
    const int recursive_adapt = 1;
    forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, do_ghost, &maxlevel);
  }

  void
  TearDown () override
  {
    t8_forest_unref (&forest);
    t8_forest_unref (&forest_adapt);
  }

  t8_forest_t forest { NULL };       /**< Uniform forest. */
  t8_forest_t forest_adapt { NULL }; /**< Forest that is adapted from \a forest. */
};

/** Check t8_forest_element_is_leaf and t8_forest_element_is_leaf_or_ghost for all leaf elements of a forest.
 * \param [in] forest  A committed forest.
 */
static void
t8_test_element_is_leaf_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Allocate memory to build a non-leaf element. */
    t8_element_t *not_leaf;
    scheme->element_new (tree_class, 1, &not_leaf);
    /* Iterate over all the tree's leaf elements, check whether the leaf
     * is correctly identified by t8_forest_element_is_leaf and t8_forest_element_is_leaf_or_ghost,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_leaf and t8_forest_element_is_leaf_or_ghost returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *leaf_element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      EXPECT_TRUE (t8_forest_element_is_leaf (forest, leaf_element, itree));
      EXPECT_TRUE (t8_forest_element_is_leaf_or_ghost (forest, leaf_element, itree, 0));
      /* Compute parent and first child of element and check that they are not in the tree */
      const int element_level = scheme->element_get_level (tree_class, leaf_element);
      if (element_level > 0) {
        scheme->element_get_parent (tree_class, leaf_element, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
        EXPECT_FALSE (t8_forest_element_is_leaf_or_ghost (forest, not_leaf, itree, 0));
      }
      if (scheme->element_is_refinable (tree_class, leaf_element)) {
        scheme->element_get_child (tree_class, leaf_element, 0, not_leaf);
        EXPECT_FALSE (t8_forest_element_is_leaf (forest, not_leaf, itree));
        EXPECT_FALSE (t8_forest_element_is_leaf_or_ghost (forest, not_leaf, itree, 0));
      }
    }
    scheme->element_destroy (tree_class, 1, &not_leaf);
  }
}

TEST_P (element_is_leaf, uniform_and_adapted_forest)
{
  t8_test_element_is_leaf_for_forest (forest);
  t8_test_element_is_leaf_for_forest (forest_adapt);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_leaf, element_is_leaf,
                          testing::Combine (AllSchemes, testing::Range (0, T8_IS_LEAF_MAX_LVL)),
                          pretty_print_eclass_scheme_and_level);

struct element_is_ghost: public testing::TestWithParam<std::tuple<int, int, cmesh_example_base *>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    const int level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = std::get<2> (GetParam ())->cmesh_create ();
    if (t8_cmesh_is_empty (cmesh)) {
      /* forest_commit does not support empty cmeshes, we skip this case */
      t8_cmesh_unref (&cmesh);
      GTEST_SKIP ();
    }
    const t8_scheme *scheme = create_from_scheme_id (scheme_id);

    const int do_ghost = 1;
    forest = t8_forest_new_uniform (cmesh, scheme, level, do_ghost, sc_MPI_COMM_WORLD);
    t8_forest_ref (forest);
    int maxlevel = level + 2;
    const int recursive_adapt = 1;
    forest_adapt = t8_forest_new_adapt (forest, t8_test_adapt_first_child, recursive_adapt, do_ghost, &maxlevel);
  }

  void
  TearDown () override
  {
    if (forest != NULL) {
      t8_forest_unref (&forest);
    }
    if (forest_adapt != NULL) {
      t8_forest_unref (&forest_adapt);
    }
  }

  t8_forest_t forest { NULL };       /**< Uniform forest with ghost layer. */
  t8_forest_t forest_adapt { NULL }; /**< Forest with ghost layer that is adapted from \a forest. */
};

/** Check t8_forest_element_is_ghost and t8_forest_element_is_leaf_or_ghost for all ghost elements of a forest.
 * \param [in] forest  A committed forest with a ghost layer.
 */
static void
t8_test_element_is_ghost_for_forest (t8_forest_t forest)
{
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  for (t8_locidx_t ighost_tree = 0; ighost_tree < num_ghost_trees; ++ighost_tree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_ghost_tree_num_leaf_elements (forest, ighost_tree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ighost_tree + num_local_trees);
    /* Allocate memory to build a non-ghost element. */
    t8_element_t *not_ghost;
    scheme->element_new (tree_class, 1, &not_ghost);
    /* Iterate over all the tree's ghost elements, check whether the ghost
     * is correctly identified by t8_forest_element_is_ghost and t8_forest_element_is_leaf_or_ghost,
     * build its parent and its first child (if they exist), and verify
     * that t8_forest_element_is_ghost and t8_forest_element_is_leaf_or_ghost returns false. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *ghost_element = t8_forest_ghost_get_leaf_element (forest, ighost_tree, ielement);
      EXPECT_TRUE (t8_forest_element_is_ghost (forest, ghost_element, ighost_tree));
      EXPECT_TRUE (t8_forest_element_is_leaf_or_ghost (forest, ghost_element, ighost_tree, 1));
      /* Compute parent and first child of element and check that they are not in the tree */
      const int element_level = scheme->element_get_level (tree_class, ghost_element);
      if (element_level > 0) {
        scheme->element_get_parent (tree_class, ghost_element, not_ghost);
        EXPECT_FALSE (t8_forest_element_is_ghost (forest, not_ghost, ighost_tree));
        EXPECT_FALSE (t8_forest_element_is_leaf_or_ghost (forest, not_ghost, ighost_tree, 1));
      }
      if (scheme->element_is_refinable (tree_class, ghost_element)) {
        scheme->element_get_child (tree_class, ghost_element, 0, not_ghost);
        EXPECT_FALSE (t8_forest_element_is_ghost (forest, not_ghost, ighost_tree));
        EXPECT_FALSE (t8_forest_element_is_leaf_or_ghost (forest, not_ghost, ighost_tree, 1));
      }
    }
    scheme->element_destroy (tree_class, 1, &not_ghost);
  }
}

TEST_P (element_is_ghost, uniform_and_adapted_forest)
{
  t8_test_element_is_ghost_for_forest (forest);
  t8_test_element_is_ghost_for_forest (forest_adapt);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_is_ghost, element_is_ghost,
                          testing::Combine (AllSchemeCollections, testing::Range (0, T8_IS_LEAF_MAX_LVL),
                                            AllCmeshsParam),
                          pretty_print_level_and_cmesh_params);
