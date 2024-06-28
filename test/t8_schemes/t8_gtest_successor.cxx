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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>

class class_successor: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();

    scheme = t8_scheme_new_default_cxx ();
    ts = scheme->eclass_schemes[eclass];

    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &successor);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &last);

    ts->t8_element_root (element);
    if (eclass == T8_ECLASS_VERTEX)
      GTEST_SKIP ();
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &successor);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &last);

    t8_scheme_cxx_unref (&scheme);
  }
  t8_eclass_t eclass;
  t8_eclass_scheme_c *ts;
  t8_scheme_cxx *scheme;
  t8_element_t *element;
  t8_element_t *successor;
  t8_element_t *child;
  t8_element_t *last;
};

/* Check the computation of the successor recursively. Iterate through the elements
 * via DFS. On the given maxlvl-1 the children are computeted iteratively. For
 * each child, the successor is checked.
 */
static void
t8_recursive_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, t8_element_t *last,
                        t8_eclass_scheme_c *ts, const int maxlvl)
{
  const int level = ts->t8_element_level (element);
  ASSERT_TRUE (ts->t8_element_level (element) <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  const int num_children = ts->t8_element_num_children (element);
  if (level == maxlvl - 1) {
    /* Check, if the successor of the last recursion is the first child of
     * of this element.
     */
    ts->t8_element_child (element, 0, child);
    EXPECT_ELEM_EQ (ts, child, successor);
    /*Check if the successor in this element is computed correctly */
    for (int ichild = 1; ichild < num_children; ichild++) {
      EXPECT_EQ (ts->t8_element_level (child), maxlvl);
      ts->t8_element_successor (child, successor, ts->t8_element_level (child), 0);
      ts->t8_element_child (element, ichild, child);
      EXPECT_ELEM_EQ (ts, child, successor);
    }
    /*If the iterator is the last element, the test can finish */
    if (ts->t8_element_equal (last, child)) {
      return;
    }
    /*Compute the next successor / "jump" out of the current element */
    else {
      EXPECT_EQ (ts->t8_element_level (child), maxlvl);
      ts->t8_element_successor (child, successor, ts->t8_element_level (child), 0);
    }
  }
  else {
    /*DFS run through the elements */
    for (int ichild = 0; ichild < num_children; ichild++) {
      ts->t8_element_child (element, ichild, child);
      t8_recursive_successor (child, successor, element, last, ts, maxlvl);
      ts->t8_element_parent (child, element);
    }
  }
}

/* Check the computation of the successor at the maximum level of the element.
 * Given the element of level maxlevel-2  with linear id 0 all children of
 * maximum level are computed. The successor runs through all these children.
 */
static void
t8_deep_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, t8_eclass_scheme_c *ts)
{
  int maxlvl = ts->t8_element_maxlevel ();
  int num_children = ts->t8_element_num_children (element);

  for (int ichild = 0; ichild < num_children; ichild++) {
    ts->t8_element_child (element, ichild, child);
    /* Go to the children at maximum level. */
    const int num_children_child = ts->t8_element_num_children (child);
    for (int jchild = 0; jchild < num_children_child; jchild++) {
      ts->t8_element_child (child, jchild, element);
      /* Check the computation of the successor. */
      ASSERT_TRUE (ts->t8_element_equal (element, successor)) << "Wrong Successor at Maxlvl.\n";
      /* Compute the next successor. */
      EXPECT_EQ (ts->t8_element_level (successor), maxlvl);
      ts->t8_element_successor (successor, successor, ts->t8_element_level (successor), 0);
    }
    ts->t8_element_parent (child, element);
  }
}

TEST_P (class_successor, test_recursive_and_deep_successor)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  /* Test at lower level. */
  for (int ilevel = 1; ilevel <= maxlvl; ilevel++) {
    ts->t8_element_set_linear_id (successor, ilevel, 0);
    ts->t8_element_last_descendant (element, last, ilevel);
    t8_recursive_successor (element, successor, child, last, ts, ilevel);
  }
  /* Test at Maxlevel. */
  ts->t8_element_set_linear_id (element, ts->t8_element_maxlevel () - 2, 0);
  ts->t8_element_set_linear_id (successor, ts->t8_element_maxlevel (), 0);
  t8_deep_successor (element, successor, last, ts);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_successor, class_successor, AllEclasses, print_eclass);
