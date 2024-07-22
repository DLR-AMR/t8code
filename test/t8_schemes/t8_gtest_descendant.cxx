/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
#include <t8_schemes/t8_standalone/t8_standalone_cxx.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>

/* This program tests the descendant function of an element. */

/* *INDENT-OFF* */
class class_schemes_descendant: public testing::TestWithParam<std::tuple<t8_eclass, int> > {
 protected:
  void
  SetUp () override
  {
    auto params = GetParam ();
    eclass = std::get<0> (params);
    int scheme_param = std::get<1> (params);
    switch (scheme_param) {
    case 0:
      scheme = t8_scheme_new_default_cxx ();
      break;
    case 1:
      scheme = t8_scheme_new_standalone_cxx ();
      break;
    default:
      SC_ABORT ("wrong scheme parameter!\n");
    }
    t8_debugf ("Creating scheme for eclass %i, schemepara: %i\n", eclass, scheme_param);
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &elem);
    ts->t8_element_new (1, &desc);
    ts->t8_element_new (1, &test);
    ts->t8_element_root (elem);
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &elem);
    ts->t8_element_destroy (1, &desc);
    ts->t8_element_destroy (1, &test);
    t8_scheme_cxx_unref (&scheme);
  }
#ifdef T8_ENABLE_DEBUG
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  t8_scheme_cxx *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;
  t8_element_t *elem;
  t8_element_t *desc;
  t8_element_t *test;
};

/* Test recursively if the first and last descendant of an element is
 * computed correctly. Only the descendant of elem->level + 1 is tested. 
 */
static void
t8_recursive_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_eclass_scheme_c *ts, int maxlvl)
{
  const int num_children = ts->t8_element_num_children (elem);
  const int level = ts->t8_element_level (elem);
  for (int ichild = 0; ichild < num_children; ichild++) {
    ts->t8_element_child (elem, ichild, desc);
    /* first child == first descendant. */
    if (ichild == 0) {
      ts->t8_element_first_descendant (elem, test, level + 1);
      EXPECT_ELEM_EQ (ts, desc, test);
    }
    /* last child == last descendant. */
    else if (ichild == num_children - 1) {
      ts->t8_element_last_descendant (elem, test, level + 1);
      EXPECT_ELEM_EQ (ts, desc, test);
    }
    else if (level > maxlvl) {
      t8_recursive_descendant (desc, elem, test, ts, maxlvl);
      ts->t8_element_parent (desc, elem);
    }
  }
}

/* Test, if the first descendant of an element is computed correctly over a range
 * of levels. 
 */
static void
t8_deep_first_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_eclass_scheme_c *ts, int level)
{
  const int elem_level = ts->t8_element_level (elem);
  ts->t8_element_copy (elem, test);

  for (int ilevel = elem_level; ilevel < level; ilevel++) {
    ts->t8_element_child (test, 0, desc);
    ts->t8_element_copy (desc, test);
  }
  ts->t8_element_first_descendant (elem, test, level);
  EXPECT_ELEM_EQ (ts, desc, test);
}

/* Test, if the last descendant of an element is computed correctly over a range
 * of levels.
 */
static void
t8_deep_last_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_eclass_scheme_c *ts, int level)
{
  ts->t8_element_copy (elem, test);

  /* Compute the correct element. */
  for (int ilevel = ts->t8_element_level (elem); ilevel < level; ilevel++) {
    const int num_children = ts->t8_element_num_children (test);
    ts->t8_element_child (test, num_children - 1, desc);
    ts->t8_element_copy (desc, test);
  }
  /* Check for equality. */
  ts->t8_element_last_descendant (elem, test, level);
  EXPECT_ELEM_EQ (ts, desc, test);
}

/* Test if the first and last descendant of an element are computed correctly.
 * The level between the element and the descendant is larger or equal to one.
 */
static void
t8_large_step_descendant (t8_element_t *elem, t8_element_t *desc, t8_element_t *test, t8_eclass_scheme_c *ts,
                          int maxlvl)
{
  for (int ilevel = ts->t8_element_level (elem); ilevel < maxlvl; ilevel++) {

    const int num_children = ts->t8_element_num_children (elem);
    /* Use these functions to perform the actual test. */
    t8_deep_first_descendant (elem, desc, test, ts, maxlvl);
    t8_deep_last_descendant (elem, desc, test, ts, maxlvl);
    for (int jchild = 0; jchild < num_children; jchild++) {
      ts->t8_element_child (elem, jchild, desc);
      t8_large_step_descendant (desc, elem, test, ts, maxlvl);
      ts->t8_element_parent (desc, elem);
    }
  }
}

TEST_P (class_schemes_descendant, test_recursive_descendant)
{
  t8_recursive_descendant (elem, desc, test, ts, maxlvl);
  t8_deep_first_descendant (elem, desc, test, ts, ts->t8_element_maxlevel ());
  t8_deep_last_descendant (elem, desc, test, ts, ts->t8_element_maxlevel ());
  t8_large_step_descendant (elem, desc, test, ts, maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_descendant, class_schemes_descendant,
                          testing::Combine (testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT), testing::Range (0, 2)));
//                          testing::Combine (testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT), testing::Range (0, 2)), print_eclass);
/* *INDENT-ON* */
