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
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>

/* This program tests the descendant function of an element. */

class class_schemes_descendant: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    scheme_id = GetParam ();

    scheme = t8_scheme_all_schemes ();
    scheme->element_new (static_cast<t8_eclass_t>(scheme_id), 1, &elem);
    scheme->element_new (static_cast<t8_eclass_t>(scheme_id), 1, &desc);
    scheme->element_new (static_cast<t8_eclass_t>(scheme_id), 1, &test);
    scheme->get_root (static_cast<t8_eclass_t>(scheme_id), elem);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (static_cast<t8_eclass_t>(scheme_id), 1, &elem);
    scheme->element_destroy (static_cast<t8_eclass_t>(scheme_id), 1, &desc);
    scheme->element_destroy (static_cast<t8_eclass_t>(scheme_id), 1, &test);
    scheme->unref ();
  }
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  int scheme_id;
  t8_scheme *scheme;
  t8_element_t *elem;
  t8_element_t *desc;
  t8_element_t *test;
};

/* Test recursively if the first and last descendant of an element is
 * computed correctly. Only the descendant of elem->level + 1 is tested. 
 */
static void
t8_recursive_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                         const int scheme_id, const int maxlvl)
{
  const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t>(scheme_id), elem);
  const int level = scheme->element_get_level (static_cast<t8_eclass_t>(scheme_id), elem);
  for (int ichild = 0; ichild < num_children; ichild++) {
    scheme->element_get_child (static_cast<t8_eclass_t>(scheme_id), elem, ichild, desc);
    /* first child == first descendant. */
    if (ichild == 0) {
      scheme->element_get_first_descendant (static_cast<t8_eclass_t>(scheme_id), elem, test, level + 1);
      EXPECT_ELEM_EQ (scheme, scheme_id, desc, test);
    }
    /* last child == last descendant. */
    else if (ichild == num_children - 1) {
      scheme->element_get_last_descendant (static_cast<t8_eclass_t>(scheme_id), elem, test, level + 1);
      EXPECT_ELEM_EQ (scheme, scheme_id, desc, test);
    }
    else if (level > maxlvl) {
      t8_recursive_descendant (desc, elem, test, scheme, static_cast<t8_eclass_t>(scheme_id), maxlvl);
      scheme->element_get_parent (static_cast<t8_eclass_t>(scheme_id), desc, elem);
    }
  }
}

/* Test, if the first descendant of an element is computed correctly over a range
 * of levels. 
 */
static void
t8_deep_first_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                          const int scheme_id, const int level)
{
  const int elem_level = scheme->element_get_level (static_cast<t8_eclass_t>(scheme_id), elem);
  scheme->element_copy (static_cast<t8_eclass_t>(scheme_id), elem, test);

  for (int ilevel = elem_level; ilevel < level; ilevel++) {
    scheme->element_get_child (static_cast<t8_eclass_t>(scheme_id), test, 0, desc);
    scheme->element_copy (static_cast<t8_eclass_t>(scheme_id), desc, test);
  }
  scheme->element_get_first_descendant (static_cast<t8_eclass_t>(scheme_id), elem, test, level);
  EXPECT_ELEM_EQ (scheme, scheme_id, desc, test);
}

/* Test, if the last descendant of an element is computed correctly over a range
 * of levels.
 */
static void
t8_deep_last_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                         const int scheme_id, const int level)
{
  scheme->element_copy (static_cast<t8_eclass_t>(scheme_id), elem, test);

  /* Compute the correct element. */
  for (int ilevel = scheme->element_get_level (static_cast<t8_eclass_t>(scheme_id), elem); ilevel < level; ilevel++) {
    const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t>(scheme_id), test);
    scheme->element_get_child (static_cast<t8_eclass_t>(scheme_id), test, num_children - 1, desc);
    scheme->element_copy (static_cast<t8_eclass_t>(scheme_id), desc, test);
  }
  /* Check for equality. */
  scheme->element_get_last_descendant (static_cast<t8_eclass_t>(scheme_id), elem, test, level);
  EXPECT_ELEM_EQ (scheme, scheme_id, desc, test);
}

/* Test if the first and last descendant of an element are computed correctly.
 * The level between the element and the descendant is larger or equal to one.
 */
static void
t8_large_step_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_scheme *scheme,
                          const int scheme_id, const int maxlvl)
{
  for (int ilevel = scheme->element_get_level (static_cast<t8_eclass_t>(scheme_id), elem); ilevel < maxlvl; ilevel++) {

    const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t>(scheme_id), elem);
    /* Use these functions to perform the actual test. */
    t8_deep_first_descendant (elem, desc, test, scheme, static_cast<t8_eclass_t>(scheme_id), maxlvl);
    t8_deep_last_descendant (elem, desc, test, scheme, static_cast<t8_eclass_t>(scheme_id), maxlvl);
    for (int jchild = 0; jchild < num_children; jchild++) {
      scheme->element_get_child (static_cast<t8_eclass_t>(scheme_id), elem, jchild, desc);
      t8_large_step_descendant (desc, elem, test, scheme, static_cast<t8_eclass_t>(scheme_id), maxlvl);
      scheme->element_get_parent (static_cast<t8_eclass_t>(scheme_id), desc, elem);
    }
  }
}

TEST_P (class_schemes_descendant, test_recursive_descendant)
{
  const int elem_maxlvl = scheme->get_maxlevel (static_cast<t8_eclass_t>(scheme_id));
  t8_recursive_descendant (elem, desc, test, scheme, scheme_id, maxlvl);
  t8_deep_first_descendant (elem, desc, test, scheme, scheme_id, elem_maxlvl);
  t8_deep_last_descendant (elem, desc, test, scheme, scheme_id, elem_maxlvl);
  t8_large_step_descendant (elem, desc, test, scheme, scheme_id, maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_descendant, class_schemes_descendant, AllSchemes);
