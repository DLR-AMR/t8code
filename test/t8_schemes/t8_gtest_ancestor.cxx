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
    scheme->element_new (eclass, 1, &correct_ancestor);
    scheme->element_new (eclass, 1, &desc_a);
    scheme->element_new (eclass, 1, &check);
    scheme->element_set_linear_id (eclass, correct_ancestor, 0, 0);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &correct_ancestor);
    scheme->element_destroy (eclass, 1, &desc_a);
    scheme->element_destroy (eclass, 1, &check);
    scheme->unref ();
  }
  /* correct_nca  -> the nearest common ancestor that we check for
    * desc_a       -> a descendant of correct_nca
    * desc_b       -> another descendant of correct_nca, different from desc_a
    * check        -> the computed nca of desc_a and desc_b, should be equal to correct_nca
    */
  t8_element_t *correct_ancestor, *desc_a, *check;
  t8_scheme *scheme;
  t8_eclass_t eclass;
};

/*Test root and parent*/
static void
t8_recursive_ancestor (t8_element_t *element, t8_element_t *child, t8_element_t *parent, t8_element_t *test_ancestor,
                       const t8_scheme *scheme, const t8_eclass_t tree_class, const int maxlvl)
{
  int num_children, i;
  const int level = scheme->element_get_level (tree_class, parent);
  const int elem_lvl = scheme->element_get_level (tree_class, element);
  T8_ASSERT (level >= elem_lvl);
  num_children = scheme->element_get_num_children (tree_class, parent);
  if (level == maxlvl) {
    return;
  }
  for (i = 0; i < num_children; i++) {
    scheme->element_get_child (tree_class, parent, i, child);
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, level, (t8_dpyramid_t *) test_ancestor);
    EXPECT_ELEM_EQ (scheme, tree_class, parent, test_ancestor);
    t8_dpyramid_ancestor ((t8_dpyramid_t *) child, elem_lvl, (t8_dpyramid_t *) test_ancestor);
    EXPECT_ELEM_EQ (scheme, tree_class, element, test_ancestor);
    t8_recursive_ancestor (element, parent, child, test_ancestor, scheme, tree_class, maxlvl);
    scheme->element_get_parent (tree_class, child, parent);
  }
}

TEST_P (ancestor, root_recursive_check)
{
  t8_element_t *parent;
  const int max_lvl = 5;
  scheme->element_new (eclass, 1, &parent);
  scheme->element_copy (eclass, correct_ancestor, parent);
  t8_recursive_ancestor (correct_ancestor, desc_a, parent, check, scheme, eclass, max_lvl);
  scheme->element_destroy (eclass, 1, &parent);
}

TEST_P (ancestor, multi_level_recursive_check)
{
  t8_element_t *parent;
  t8_element_t *correct_ancestor_high_level;
  int recursion_depth = 5;
  const int max_lvl = scheme->get_maxlevel (eclass);
  int i;
  scheme->element_new (eclass, 1, &parent);
  scheme->element_new (eclass, 1, &correct_ancestor_high_level);

  t8_gloidx_t leaves_on_level;
  for (i = recursion_depth; i < max_lvl; i++) {
    leaves_on_level = scheme->element_count_leaves (eclass, correct_ancestor, i - recursion_depth);
    scheme->element_set_linear_id (eclass, correct_ancestor_high_level, i - recursion_depth, leaves_on_level / 2);
    scheme->element_copy (eclass, correct_ancestor_high_level, parent);
    t8_recursive_ancestor (correct_ancestor, desc_a, parent, check, scheme, eclass, i);
  }
  scheme->element_destroy (eclass, 1, &parent);
  scheme->element_destroy (eclass, 1, &correct_ancestor_high_level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ancestorestor, ancestor, testing::Values (T8_ECLASS_PYRAMID), print_eclass);
