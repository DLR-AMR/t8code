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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>

/* This program tests the descendant function of an element. */

class class_schemes_descendant: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    tree_class = GetParam ();

    scheme = t8_scheme_new_default ();
    scheme->element_new (tree_class, 1, &elem);
    scheme->element_new (tree_class, 1, &desc);
    scheme->element_new (tree_class, 1, &test);
    scheme->get_root (tree_class, elem);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (tree_class, 1, &elem);
    scheme->element_destroy (tree_class, 1, &desc);
    scheme->element_destroy (tree_class, 1, &test);
    scheme->unref ();
  }
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  t8_scheme *scheme;
  t8_eclass_t tree_class;
  t8_element_t *elem;
  t8_element_t *desc;
  t8_element_t *test;
};

/* Test recursively if the first and last descendant of an element is
 * computed correctly. Only the descendant of elem->level + 1 is tested. 
 */
static void
t8_recursive_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                         const t8_eclass_t tree_class, int maxlvl)
{
  const int num_children = scheme->element_get_num_children (tree_class, elem);
  const int level = scheme->element_get_level (tree_class, elem);
  for (int ichild = 0; ichild < num_children; ichild++) {
    scheme->element_get_child (tree_class, elem, ichild, desc);
    /* first child == first descendant. */
    if (ichild == 0) {
      scheme->element_get_first_descendant (tree_class, elem, test, level + 1);
      EXPECT_ELEM_EQ (scheme, tree_class, desc, test);
    }
    /* last child == last descendant. */
    else if (ichild == num_children - 1) {
      scheme->element_get_last_descendant (tree_class, elem, test, level + 1);
      EXPECT_ELEM_EQ (scheme, tree_class, desc, test);
    }
    else if (level > maxlvl) {
      t8_recursive_descendant (desc, elem, test, scheme, tree_class, maxlvl);
      scheme->element_get_parent (tree_class, desc, elem);
    }
  }
}

/* Test, if the first descendant of an element is computed correctly over a range
 * of levels. 
 */
static void
t8_deep_first_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                          const t8_eclass_t tree_class, int level)
{
  const int elem_level = scheme->element_get_level (tree_class, elem);
  scheme->element_copy (tree_class, elem, test);

  for (int ilevel = elem_level; ilevel < level; ilevel++) {
    scheme->element_get_child (tree_class, test, 0, desc);
    scheme->element_copy (tree_class, desc, test);
  }
  scheme->element_get_first_descendant (tree_class, elem, test, level);
  EXPECT_ELEM_EQ (scheme, tree_class, desc, test);
}

/* Test, if the last descendant of an element is computed correctly over a range
 * of levels.
 */
static void
t8_deep_last_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                         const t8_eclass_t tree_class, int level)
{
  scheme->element_copy (tree_class, elem, test);

  /* Compute the correct element. */
  for (int ilevel = scheme->element_get_level (tree_class, elem); ilevel < level; ilevel++) {
    const int num_children = scheme->element_get_num_children (tree_class, test);
    scheme->element_get_child (tree_class, test, num_children - 1, desc);
    scheme->element_copy (tree_class, desc, test);
  }
  /* Check for equality. */
  scheme->element_get_last_descendant (tree_class, elem, test, level);
  EXPECT_ELEM_EQ (scheme, tree_class, desc, test);
}

/* Test if the first and last descendant of an element are computed correctly.
 * The level between the element and the descendant is larger or equal to one.
 */
static void
t8_large_step_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                          const t8_eclass_t tree_class, int maxlvl)
{
  for (int ilevel = scheme->element_get_level (tree_class, elem); ilevel < maxlvl; ilevel++) {

    const int num_children = scheme->element_get_num_children (tree_class, elem);
    /* Use these functions to perform the actual test. */
    t8_deep_first_descendant (elem, desc, test, scheme, tree_class, maxlvl);
    t8_deep_last_descendant (elem, desc, test, scheme, tree_class, maxlvl);
    for (int jchild = 0; jchild < num_children; jchild++) {
      scheme->element_get_child (tree_class, elem, jchild, desc);
      t8_large_step_descendant (desc, elem, test, scheme, tree_class, maxlvl);
      scheme->element_get_parent (tree_class, desc, elem);
    }
  }
}

TEST_P (class_schemes_descendant, test_recursive_descendant)
{
  t8_recursive_descendant (elem, desc, test, scheme, tree_class, maxlvl);
  t8_deep_first_descendant (elem, desc, test, scheme, tree_class, scheme->get_maxlevel (tree_class));
  t8_deep_last_descendant (elem, desc, test, scheme, tree_class, scheme->get_maxlevel (tree_class));
  t8_large_step_descendant (elem, desc, test, scheme, tree_class, maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_descendant, class_schemes_descendant, AllEclasses, print_eclass);
