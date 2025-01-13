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
#include "t8_gtest_dfs_base.hxx"
#include <test/t8_gtest_macros.hxx>

class class_find_parent: public TestDFS {
  void
  check_element () override
  {
    const int num_children = scheme->element_get_num_children (eclass, element);
    for (int ichild = 0; ichild < num_children; ichild++) {
      scheme->element_get_child (eclass, element, ichild, child);
      /* Compute parent of child */
      scheme->element_get_parent (eclass, child, test_parent);
      /* Check that it is equal to the original element */
      EXPECT_ELEM_EQ (scheme, eclass, element, test_parent);
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    scheme->element_new (eclass, 1, &child);
    scheme->element_new (eclass, 1, &test_parent);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    scheme->element_destroy (eclass, 1, &child);
    scheme->element_destroy (eclass, 1, &test_parent);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t *child;
  t8_element_t *test_parent;
};

TEST_P (class_find_parent, t8_compute_child_find_parent)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 4;
#else
  const int maxlvl = 6;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_find_parent, class_find_parent, AllSchemes);
