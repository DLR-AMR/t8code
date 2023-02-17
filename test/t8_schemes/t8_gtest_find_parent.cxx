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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* *INDENT-OFF* */
class class_find_parent:public testing::TestWithParam <t8_eclass_t > {
protected:
  void SetUp () override {
    eclass = GetParam();

    scheme = t8_scheme_new_default_cxx ();
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
  }
  void TearDown () override {
    /* Destroy scheme */
    t8_scheme_cxx_unref (&scheme);
  }
  t8_eclass_t           eclass;
  t8_scheme_cxx         *scheme;
  t8_eclass_scheme_c    *ts;

};
/* *INDENT-ON* */

static void
t8_recursive_child_find_parent (t8_element_t *element, t8_element_t *child,
                                t8_element_t *test_parent,
                                t8_eclass_scheme_c *ts, int level,
                                const int maxlvl)
{
  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);

  /* Get number of children */
  int                 num_children = ts->t8_element_num_children (element);
  /* Get child and test_parent, to check if test_parent = parent of child */
  if (level == maxlvl)
    return;
  for (int ichild = 0; ichild < num_children; ichild++) {
    /* Compute child i */
    ts->t8_element_child (element, ichild, child);
    /* Compute parent of child */
    ts->t8_element_parent (child, test_parent);
    /* If its equal, call child_find_parent, to check if parent-child relation
     * is correct in next level until maxlvl is reached*/
    ASSERT_TRUE (!ts->
                 t8_element_compare (element,
                                     test_parent)) <<
      "Computed child_parent is not the parent.";

    t8_recursive_child_find_parent (child, element, test_parent, ts,
                                    level + 1, maxlvl);
    /* After the check we know the parent-function is correct for this child.
     * Therefore we can use it to recompute the element*/
    ts->t8_element_parent (child, element);
  }
}

TEST_P (class_find_parent, t8_compute_child_find_parent)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int           maxlvl = 4;
#else
  const int           maxlvl = 6;
#endif

  t8_element_t       *element;
  t8_element_t       *child;
  t8_element_t       *test_parent;

  /* Get element and initialize it */
  ts->t8_element_new (1, &element);
  ts->t8_element_new (1, &child);
  ts->t8_element_new (1, &test_parent);

  ts->t8_element_set_linear_id (element, 0, 0);
  /* Check for correct parent-child relation */
  t8_recursive_child_find_parent (element, child, test_parent, ts, 0, maxlvl);
  /* Destroy element */
  ts->t8_element_destroy (1, &element);
  ts->t8_element_destroy (1, &child);
  ts->t8_element_destroy (1, &test_parent);
}

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_find_parent, class_find_parent,testing::Range(T8_ECLASS_ZERO, T8_ECLASS_COUNT));
/* *INDENT-ON* */
