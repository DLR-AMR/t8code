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

/** \file t8_gtest_nca.cxx
* Provide tests to check the functionality of the nearest-common-ancestor function
* for every element.
*/

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

class nca: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    tree_class = GetParam ();
    scheme = t8_scheme_new_default ();
    scheme->element_new (tree_class, 1, &correct_nca);
    scheme->element_new (tree_class, 1, &desc_a);
    scheme->element_new (tree_class, 1, &desc_b);
    scheme->element_new (tree_class, 1, &check);
    scheme->element_set_linear_id (tree_class, correct_nca, 0, 0);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (tree_class, 1, &correct_nca);
    scheme->element_destroy (tree_class, 1, &desc_a);
    scheme->element_destroy (tree_class, 1, &desc_b);
    scheme->element_destroy (tree_class, 1, &check);
    scheme->unref ();
  }
  /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
  t8_element_t *correct_nca, *desc_a, *desc_b, *check;
  t8_scheme *scheme;
  t8_eclass_t tree_class;
};

/**
 * Test the nca for the children of the root-element
 * 
 */
TEST_P (nca, nca_check_shallow)
{
  int i, j;
  const int num_children = scheme->element_get_num_children (tree_class, correct_nca);
  /* Iterate over all combinations of two children from correct_nca */
  for (i = 0; i < num_children - 1; i++) {
    scheme->element_get_child (tree_class, correct_nca, i, desc_a);
    for (j = i + 1; j < num_children; j++) {
      scheme->element_get_child (tree_class, correct_nca, j, desc_b);
      /*Compute the nca */
      scheme->element_get_nca (tree_class, desc_a, desc_b, check);
      /*expect equality */
      EXPECT_ELEM_EQ (scheme, tree_class, check, correct_nca);
    }
  }
}

/**
 * Check the nca for elements on the maximal level.
 * We iteratively compute an element that is the correct nca up to level 10.
 * The first and the last descendant of this element at the maximal level are
 * computed and used as an input for element_get_nca.
 */
TEST_P (nca, nca_check_deep)
{
  const int max_lvl = 10;
  const int elem_max_level = scheme->get_maxlevel (tree_class);
  /* num_children is not a const here, cause this can change for pyramids */
  int num_children;
  /* iterate over levels and children */
  int lvl, check_lvl_a, check_lvl_b, child_id;
  t8_element_t *tmp;

  scheme->element_new (tree_class, 1, &tmp);
  scheme->element_copy (tree_class, correct_nca, tmp);
  for (lvl = 1; lvl <= max_lvl; lvl++) {
    num_children = scheme->element_get_num_children (tree_class, tmp);
    for (child_id = 0; child_id < num_children; child_id++) {
      scheme->element_get_child (tree_class, tmp, child_id, correct_nca);
      /* Compute first and last descendant at every level up to elem_max_lvl. 
       * They have the correct_nca as the nca */
      for (check_lvl_a = lvl + 1; check_lvl_a < elem_max_level; check_lvl_a++) {
        scheme->element_construct_first_descendant (tree_class, correct_nca, desc_a, check_lvl_a);
        for (check_lvl_b = lvl + 1; check_lvl_b < elem_max_level; check_lvl_b++) {
          scheme->element_construct_last_descendant (tree_class, correct_nca, desc_b, check_lvl_b);
          /* Compute the nca of desc_a and desc_b */
          scheme->element_get_nca (tree_class, desc_a, desc_b, check);
          if (tree_class == T8_ECLASS_VERTEX) {
            /* first- last-descendant logic does not hold for vertices. */
            EXPECT_EQ (
              scheme->element_get_level (tree_class, check),
              SC_MIN (scheme->element_get_level (tree_class, desc_a), scheme->element_get_level (tree_class, desc_b)));
          }
          else {
            /* Expect equality of correct_nca and check for every other class */
            EXPECT_ELEM_EQ (scheme, tree_class, correct_nca, check);
          }
        }
      }
    }
    /* Determine next element */
    if (num_children != 1) {
      scheme->element_copy (tree_class, tmp, correct_nca);
      /* Continue in the middle */
      scheme->element_get_child (tree_class, correct_nca, num_children / 2, tmp);
    }
    else {
      scheme->element_copy (tree_class, correct_nca, tmp);
    }
  }
  scheme->element_destroy (tree_class, 1, &tmp);
}

/**
 * Recursively check the computation of the nca of all possible combination of descendants of the
 * \a correct_nca that have \a correct_nca as the nca. 
 * 
 * \param[in] correct_nca       The correct nearest common ancestor
 * \param[in] desc_a            Storage for the computation of a descendant of \correct_nca
 * \param[in] desc_b            Storage for the computation of a descendant of \correct_nca
 * \param[in] check             Storage for the computation of the nca of \a desc_a and \a desc_b
 * \param[in] parent_a          An initialized element, descendant of \a correct_nca, not a descendant or ancestor of \a parent_b. \a desc_a will be a child of it
 * \param[in] parent_b          An initialized element, descendant of \a correct_nca, not a descendant or ancestor of \a parent_a. \a desc_b will be a child of it
 * \param[in] max_lvl           The maximal depth of the recursion
 * \param[in] ts                the scheme to use
 */
static void
t8_recursive_nca_check (t8_element_t *check_nca, t8_element_t *desc_a, t8_element_t *desc_b, t8_element_t *check,
                        t8_element_t *parent_a, t8_element_t *parent_b, const int max_lvl, t8_scheme *scheme,
                        const t8_eclass_t tree_class)
{
  T8_ASSERT (max_lvl <= scheme->get_maxlevel (tree_class) - 1);
  /* compute the level of the parents */
  int level_a = scheme->element_get_level (tree_class, parent_a);
  int level_b = scheme->element_get_level (tree_class, parent_b);
  int num_children_a, num_children_b;
  int i, j;
  /* If one parent has reached the maximal level, the test returns */
  if (level_a == max_lvl || level_b == max_lvl) {
    return;
  }
  num_children_a = scheme->element_get_num_children (tree_class, parent_a);
  num_children_b = scheme->element_get_num_children (tree_class, parent_b);

  /* Iterate over all children of parent_a */
  for (i = 0; i < num_children_a; i++) {
    scheme->element_get_child (tree_class, parent_a, i, desc_a);
    /* Iterate over all children of parent_b */
    for (j = 0; j < num_children_b; j++) {
      scheme->element_get_child (tree_class, parent_b, j, desc_b);
      scheme->element_get_nca (tree_class, desc_a, desc_b, check);

      if (!scheme->element_is_equal (tree_class, check_nca, check)) {
        level_a = scheme->element_get_level (tree_class, desc_a);
        level_b = scheme->element_get_level (tree_class, desc_b);

        int level_c = scheme->element_get_level (tree_class, check_nca);
        int level_nca = scheme->element_get_level (tree_class, check);
        /* Output the linear id of the descendants where the computation fails.
         * This makes debugging a lot easier, as one can reconstruct the descendants
         * via t8_element_set_linear_id and can directly test them instead of waiting
         * until the recursion reaches the faulty computation. */
        t8_debugf ("id of desc_a: %li, level: %i\n",
                   static_cast<long> (scheme->element_get_linear_id (tree_class, desc_a, level_a)), level_a);
        t8_debugf ("id of desc_b: %li, level: %i\n",
                   static_cast<long> (scheme->element_get_linear_id (tree_class, desc_b, level_b)), level_b);

        for (int k = SC_MAX (level_a, level_b); k >= 0; k--) {
          t8_debugf ("id of desc_a: %li, level: %i\n",
                     static_cast<long> (scheme->element_get_linear_id (tree_class, desc_a, k)), k);
          t8_debugf ("id of desc_b: %li, level: %i\n",
                     static_cast<long> (scheme->element_get_linear_id (tree_class, desc_b, k)), k);
        }

        t8_debugf ("id of the correct nca: %li, level: %i\n",
                   static_cast<long> (scheme->element_get_linear_id (tree_class, check_nca, level_c)), level_c);

        t8_debugf ("id of the computed nca: %li, level: %i\n",
                   static_cast<long> (scheme->element_get_linear_id (tree_class, check, level_nca)), level_nca);

        SC_ABORT ("Computed nca is not the correct nca!\n");
      }
      /* parent_a stays fixed, b-part goes one level deeper into the recursion */
      t8_recursive_nca_check (check_nca, desc_a, parent_b, check, parent_a, desc_b, max_lvl, scheme, tree_class);
      /* We reused parent_b, hence we have to recompute the correct parent */
      scheme->element_get_parent (tree_class, desc_b, parent_b);
    }
    /* We reused parent_a, hence we have to recompute the correct parent */
    scheme->element_get_parent (tree_class, desc_a, parent_a);
  }
}

/* Recursively check the computation of the nca. recursion_depth defines up to which
 * level we compute descendants of correct_nca that should have correct_nca as the 
 * output of element_get_nca.*/
TEST_P (nca, recursive_check)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int recursion_depth = 3;
#else
  /* User lower recursion depth for pyramids, it takes to much time otherwise */
  const int recursion_depth = 4;
#endif
  t8_element_t *parent_a, *parent_b;
  int num_children;
  num_children = scheme->element_get_num_children (tree_class, correct_nca);
  int i, j;
  if (num_children > 1) {
    scheme->element_new (tree_class, 1, &parent_a);
    scheme->element_new (tree_class, 1, &parent_b);
    scheme->element_get_child (tree_class, correct_nca, 0, parent_a);
    scheme->element_get_child (tree_class, correct_nca, 1, parent_b);
    for (i = 0; i < num_children - 1; i++) {
      scheme->element_get_child (tree_class, correct_nca, i, parent_a);
      for (j = i + 1; j < num_children; j++) {
        scheme->element_get_child (tree_class, correct_nca, j, parent_b);
        t8_recursive_nca_check (correct_nca, desc_a, desc_b, check, parent_a, parent_b, recursion_depth, scheme,
                                tree_class);
      }
    }
  }
  else {
    GTEST_SKIP ();
  }
  scheme->element_destroy (tree_class, 1, &parent_a);
  scheme->element_destroy (tree_class, 1, &parent_b);
}

/* Test the nca recursively for elements in the middle of the uniform refinement tree
 * up to the maximal level. 
 * Be careful when increasing the recursion_depth, as it increases the number of test-cases exponentially. */
TEST_P (nca, recursive_check_higher_level)
{

  const int recursion_depth = 3;
  const int max_lvl = scheme->get_maxlevel (tree_class);
  t8_element_t *parent_a;
  t8_element_t *parent_b;
  t8_element_t *correct_nca_high_level;
  int num_children;
  int i, k, l;
  t8_gloidx_t leaves_on_level;
  EXPECT_TRUE (max_lvl - recursion_depth >= 0);

  scheme->element_new (tree_class, 1, &parent_a);
  scheme->element_new (tree_class, 1, &parent_b);
  scheme->element_new (tree_class, 1, &correct_nca_high_level);

  /* Test on different levels around the middle of the refinement tree */
  for (i = recursion_depth; i < max_lvl; i++) {
    leaves_on_level = scheme->element_count_leaves (tree_class, correct_nca, i - recursion_depth);
    /* middle = leaves/2 */
    scheme->element_set_linear_id (tree_class, correct_nca_high_level, i - recursion_depth, leaves_on_level / 2);

    /* Initialization for recursive_nca_check */
    num_children = scheme->element_get_num_children (tree_class, correct_nca_high_level);
    if (num_children > 1) {
      /* Compute children on to different branches in the tree an test them. 
       * This ensures, that the nca of all their descendants has to be correct_nca_high_level*/
      for (k = 0; k < num_children; k++) {
        scheme->element_get_child (tree_class, correct_nca_high_level, k, parent_a);
        for (l = 0; l < num_children; l++) {
          scheme->element_get_child (tree_class, correct_nca_high_level, l, parent_b);
          if (k != l) {
            t8_recursive_nca_check (correct_nca_high_level, desc_a, desc_b, check, parent_a, parent_b, i, scheme,
                                    tree_class);
          }
          else {
            scheme->element_get_nca (tree_class, parent_a, parent_b, check);
            EXPECT_ELEM_EQ (scheme, tree_class, parent_a, check);
            EXPECT_ELEM_EQ (scheme, tree_class, parent_b, check);
          }
        }
      }
    }
    else {
      scheme->element_destroy (tree_class, 1, &parent_a);
      scheme->element_destroy (tree_class, 1, &parent_b);
      scheme->element_destroy (tree_class, 1, &correct_nca_high_level);
      GTEST_SKIP ();
    }
  }
  /* Clean-up */
  scheme->element_destroy (tree_class, 1, &parent_a);
  scheme->element_destroy (tree_class, 1, &parent_b);
  scheme->element_destroy (tree_class, 1, &correct_nca_high_level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_nca, nca, AllEclasses, print_eclass);
