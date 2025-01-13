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
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"

class class_test_equal: public TestDFS {
  void
  check_element () override
  {
    const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t> (scheme_id), element);
    for (int ichild1 = 0; ichild1 < num_children; ichild1++) {
      scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, ichild1, child1);
      /* the child must be different than its parent */
      EXPECT_FALSE (scheme->element_is_equal (static_cast<t8_eclass_t> (scheme_id), element, child1));
      for (int ichild2 = 0; ichild2 < num_children; ichild2++) {
        scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, ichild2, child2);
        /* the child must be different than its parent */
        EXPECT_FALSE (scheme->element_is_equal (static_cast<t8_eclass_t> (scheme_id), element, child2));
        const int equal = scheme->element_is_equal (static_cast<t8_eclass_t> (scheme_id), child1, child2);
        /* The children must be equal if and only if their indices are equal. */
        EXPECT_EQ (equal, ichild1 == ichild2);
        /* t8_element_equal should compute the same as t8_element_compare, 
         * when we only check if compare has 0 as result. */
        const int compare_equal = !scheme->element_compare (static_cast<t8_eclass_t> (scheme_id), child1, child2);
        EXPECT_EQ (equal, compare_equal);
      }
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &child1);
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &child2);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &child1);
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &child2);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t *child1;
  t8_element_t *child2;
};

TEST_P (class_test_equal, test_equal_dfs)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, AllSchemes);
