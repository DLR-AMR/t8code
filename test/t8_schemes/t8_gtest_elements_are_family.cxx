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

#include <gtest/gtest.h>
#include <t8_eclass/t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"

struct are_family: public TestDFS
{
 private:
  void
  check_element () override
  {
    if (eclass == T8_ECLASS_VERTEX) {
      GTEST_SKIP ();
    }
    const int num_children = scheme->element_get_num_children (eclass, element);
    t8_element_t **element_array = T8_TESTSUITE_ALLOC (t8_element_t *, num_children);
    scheme->element_new (eclass, num_children, element_array);
    scheme->element_get_children (eclass, element, num_children, element_array);
    EXPECT_TRUE (scheme->elements_are_family (eclass, element_array));
    scheme->element_get_child (eclass, element_array[num_children - 1], 0, element_array[num_children - 1]);
    EXPECT_FALSE (scheme->elements_are_family (eclass, element_array));

    scheme->element_destroy (eclass, num_children, element_array);
    T8_TESTSUITE_FREE (element_array);
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

TEST_P (are_family, test_are_family)
{
#if T8_TEST_LEVEL_INT >= 1
  const int maxlvl = 4;
#else
  const int maxlvl = 6;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_elements_are_family, are_family, AllSchemes, print_all_schemes);
