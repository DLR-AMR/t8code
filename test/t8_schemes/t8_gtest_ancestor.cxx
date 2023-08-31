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
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

class ancestor: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    scheme = t8_scheme_new_default_cxx ();
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &correct_anc);
    ts->t8_element_new (1, &desc_a);
    ts->t8_element_new (1, &check);
    ts->t8_element_set_linear_id (correct_anc, 0, 0);
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &correct_anc);
    ts->t8_element_destroy (1, &desc_a);
    ts->t8_element_destroy (1, &check);
    t8_scheme_cxx_unref (&scheme);
  }
  /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
  t8_element_t *correct_anc, *desc_a, *check;
  t8_scheme_cxx *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;
};

/*Test root and parent*/
static void
t8_recursive_ancestor (t8_element_t *element, t8_element_t *child, t8_element_t *parent, t8_element_t *test_anc,
                       t8_eclass_scheme_c *ts, const int maxlvl)
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
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, level, (t8_dpyramid_t *) test_anc);
    SC_CHECK_ABORT (!ts->t8_element_compare (parent, test_anc), "Computed ancestor is not equal to the parent\n");
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, elem_lvl, (t8_dpyramid_t *) test_anc);
    SC_CHECK_ABORT (!ts->t8_element_compare (element, test_anc),
                    "Computed ancestor is not equal to the correct ancestor\n");
    t8_recursive_ancestor (element, parent, child, test_anc, ts, maxlvl);
    ts->t8_element_parent (child, parent);
  }
}

TEST_P (ancestor, root_recursive_check)
{
  t8_element_t *parent;
  int max_lvl = 5;
  ts->t8_element_new (1, &parent);
  ts->t8_element_copy (correct_anc, parent);
  t8_recursive_ancestor (correct_anc, desc_a, parent, check, ts, max_lvl);
  ts->t8_element_destroy (1, &parent);
}

TEST_P (ancestor, multi_level_recursive_check)
{
  t8_element_t *parent;
  t8_element_t *correct_anc_high_level;
  int recursion_depth = 5;
  int max_lvl = ts->t8_element_maxlevel ();
  int i;
  ts->t8_element_new (1, &parent);
  ts->t8_element_new (1, &correct_anc_high_level);

  t8_gloidx_t leafs_on_level;
  for (i = recursion_depth; i < max_lvl; i++) {
    leafs_on_level = ts->t8_element_count_leafs (correct_anc, i - recursion_depth);
    ts->t8_element_set_linear_id (correct_anc_high_level, i - recursion_depth, leafs_on_level / 2);
    ts->t8_element_copy (correct_anc_high_level, parent);
    t8_recursive_ancestor (correct_anc, desc_a, parent, check, ts, i);
  }
  ts->t8_element_destroy (1, &parent);
  ts->t8_element_destroy (1, &correct_anc_high_level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ancestor, ancestor, testing::Values (T8_ECLASS_PYRAMID));
