/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_gtest_is_ancestor.cxx
 * This test checks the element_is_ancestor function of the scheme interface.
 * Starting from an element we build up all its ancestor up to root level and 
 * test whether the element_is_ancestor function returns true.
 * We also test sone cases where the function returns false, namely when putting
 * in an element and an ancestor in reverse order.
 * 
 * This test is modified from the ancestor_id test.
 */

#include <gtest/gtest.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include "t8_gtest_dfs_base.hxx"
#include <test/t8_gtest_macros.hxx>

class class_is_ancestor: public TestDFS {
  void
  check_element () override
  {
    t8_element_t *ancestor;
    scheme->element_new (eclass, 1, &ancestor);
    /* Get level of current element */
    const int level = scheme->element_get_level (eclass, element);

    /* Iterate over all levels above the current element and check
    if ancestor id corresponds with the child id of elem, parent, grandparent, ... */
    for (int levels_above_elem = 0; levels_above_elem < level; levels_above_elem++) {
      const int ancestor_level = level - levels_above_elem;
      /* Compute the ancestor by iteratively computing the parent */
      scheme->element_copy (eclass, element, ancestor);
      for (int level_diff = 0; level_diff < levels_above_elem; level_diff++) {
        scheme->element_get_parent (eclass, ancestor, ancestor);
      }
      // Check whether element_is_ancestor correctly identifies our candidate as an ancestor
      EXPECT_TRUE (scheme->element_is_ancestor (eclass, ancestor, element));
      // We check that element_is_ancestor properly returns false by
      // checking that element is not an ancestor of our candidate if their levels do not agree.
      if (ancestor_level != level) {
        EXPECT_FALSE (scheme->element_is_ancestor (eclass, element, ancestor));
      }
    }
    scheme->element_destroy (eclass, 1, &ancestor);
  }

 protected:
  void
  SetUp () override
  {
    /* Setup DFS test */
    dfs_test_setup ();
  }
  void
  TearDown () override
  {
    /* Destroy DFS test */
    dfs_test_teardown ();
  }
};

TEST_P (class_is_ancestor, t8_recursive_dfs_is_ancestor)
{
#if T8_TEST_LEVEL_INT >= 1  // test level medium or lower
  const int maxlvl = 4;
#else
  const int maxlvl = 6;  // test level full
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_is_ancestor, class_is_ancestor, AllSchemes, print_all_schemes);
INSTANTIATE_TEST_SUITE_P (t8_gtest_is_ancestor_multilevel, class_is_ancestor, AllMultilevelSchemes, print_all_schemes);

/** Siblings are never ancestors of each other. In the multilevel schemes this especially holds for
 * a refined element, which is its own first child, and the other children.
 * The siblings are constructed as children of the parent, since not all schemes implement element_get_sibling.
 */
struct class_is_ancestor_siblings: public TestDFS
{
  void
  check_element () override
  {
    if (scheme->element_get_level (eclass, element) == 0) {
      return;
    }
    /* Only the default quad and hex schemes implement element_get_sibling. */
    const bool check_get_sibling
      = std::get<0> (GetParam ()) == 2 && (eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX);

    t8_element_t *parent;
    t8_element_t *sibling;
    t8_element_t *sibling_by_id;
    scheme->element_new (eclass, 1, &parent);
    scheme->element_new (eclass, 1, &sibling);
    scheme->element_new (eclass, 1, &sibling_by_id);
    scheme->element_get_parent (eclass, element, parent);

    const int child_id = scheme->element_get_child_id (eclass, element);
    const int num_siblings = scheme->element_get_num_siblings (eclass, element);
    EXPECT_EQ (num_siblings, scheme->element_get_num_children (eclass, parent));
    for (int sibling_id = 0; sibling_id < num_siblings; ++sibling_id) {
      scheme->element_get_child (eclass, parent, sibling_id, sibling);
      if (check_get_sibling) {
        scheme->element_get_sibling (eclass, element, sibling_id, sibling_by_id);
        EXPECT_ELEM_EQ (scheme, eclass, sibling_by_id, sibling);
      }
      if (sibling_id == child_id) {
        EXPECT_ELEM_EQ (scheme, eclass, sibling, element);
        continue;
      }
      EXPECT_FALSE (scheme->element_is_ancestor (eclass, element, sibling));
      EXPECT_FALSE (scheme->element_is_ancestor (eclass, sibling, element));
    }
    scheme->element_destroy (eclass, 1, &sibling_by_id);
    scheme->element_destroy (eclass, 1, &sibling);
    scheme->element_destroy (eclass, 1, &parent);
  }
};

TEST_P (class_is_ancestor_siblings, t8_recursive_dfs_siblings_are_no_ancestors)
{
#if T8_TEST_LEVEL_INT >= 1  // test level medium or lower
  const int maxlvl = 4;
#else
  const int maxlvl = 6;  // test level full
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_is_ancestor_siblings, class_is_ancestor_siblings, AllMultilevelSchemes,
                          print_all_schemes);
