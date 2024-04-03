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
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"

class class_test_equal: public TestDFS {
  virtual void
  check_element ()
  {
    const int num_children = ts->t8_element_num_children (element);
    for (int ichild1 = 0; ichild1 < num_children; ichild1++) {
      ts->t8_element_child (element, ichild1, child1);
      /* the child must be different than its parent */
      EXPECT_FALSE (ts->t8_element_equal (element, child1));
      for (int ichild2 = 0; ichild2 < num_children; ichild2++) {
        ts->t8_element_child (element, ichild2, child2);
        /* the child must be different than its parent */
        EXPECT_FALSE (ts->t8_element_equal (element, child2));
        const int equal = ts->t8_element_equal (child1, child2);
        /* The children must be equal if and only if their indices are equal. */
        EXPECT_EQ (equal, ichild1 == ichild2);
        /* t8_element_equal should compute the same as t8_element_compare, 
         * when we only check if compare has 0 as result. */
        const int compare_equal = !ts->t8_element_compare (child1, child2);
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
    ts->t8_element_new (1, &child1);
    ts->t8_element_new (1, &child2);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    ts->t8_element_destroy (1, &child1);
    ts->t8_element_destroy (1, &child2);

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

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, AllEclasses, print_eclass);
