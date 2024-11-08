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
#include <test/t8_gtest_macros.hxx>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_default.hxx>

class ancestor: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    scheme = t8_scheme_new_default ();
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &correct_ancestor);
    ts->t8_element_new (1, &desc_a);
    ts->t8_element_new (1, &check);
    ts->t8_element_set_linear_id (correct_ancestor, 0, 0);
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &correct_ancestor);
    ts->t8_element_destroy (1, &desc_a);
    ts->t8_element_destroy (1, &check);
    t8_schemexx_unref (&scheme);
  }
  /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
  t8_element_t *correct_ancestor, *desc_a, *check;
  t8_scheme *scheme;
  t8_scheme *ts;
  t8_eclass_t eclass;
};

/*Test root and parent*/
static void
t8_recursive_ancestor (t8_element_t *element, t8_element_t *child, t8_element_t *parent, t8_element_t *test_ancestor,
                       t8_scheme *ts, const int maxlvl)
{
  int num_children, i;
  int level = ts->t8_element_level (parent);
  int elem_lvl = ts->t8_element_level (element);
  T8_ASSERT (level >= elem_lvl);
  num_children = ts->t8_element_num_children (parent);
  if (level == maxlvl) {
    return;
  }
  for (i = 0; i < num_children; i++) {
    ts->t8_element_child (parent, i, child);
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, level, (t8_dpyramid_t *) test_ancestor);
    EXPECT_ELEM_EQ (ts, parent, test_ancestor);
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, elem_lvl, (t8_dpyramid_t *) test_ancestor);
    EXPECT_ELEM_EQ (ts, element, test_ancestor);
    t8_recursive_ancestor (element, parent, child, test_ancestor, ts, maxlvl);
    ts->t8_element_parent (child, parent);
  }
}

TEST_P (ancestor, root_recursive_check)
{
  t8_element_t *parent;
  int max_lvl = 5;
  ts->t8_element_new (1, &parent);
  ts->t8_element_copy (correct_ancestor, parent);
  t8_recursive_ancestor (correct_ancestor, desc_a, parent, check, ts, max_lvl);
  ts->t8_element_destroy (1, &parent);
}

TEST_P (ancestor, multi_level_recursive_check)
{
  t8_element_t *parent;
  t8_element_t *correct_ancestor_high_level;
  int recursion_depth = 5;
  int max_lvl = ts->t8_element_maxlevel ();
  int i;
  ts->t8_element_new (1, &parent);
  ts->t8_element_new (1, &correct_ancestor_high_level);

  t8_gloidx_t leaves_on_level;
  for (i = recursion_depth; i < max_lvl; i++) {
    leaves_on_level = ts->t8_element_count_leaves (correct_ancestor, i - recursion_depth);
    ts->t8_element_set_linear_id (correct_ancestor_high_level, i - recursion_depth, leaves_on_level / 2);
    ts->t8_element_copy (correct_ancestor_high_level, parent);
    t8_recursive_ancestor (correct_ancestor, desc_a, parent, check, ts, i);
  }
  ts->t8_element_destroy (1, &parent);
  ts->t8_element_destroy (1, &correct_ancestor_high_level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ancestorestor, ancestor, testing::Values (T8_ECLASS_PYRAMID), print_eclass);
