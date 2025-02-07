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

/** In this test we iterate through all elements. 
 * For every leaf we check if the element is equal to the element we get when setting it from the linear id.
 * The id_counter is then increased to match the id of the next leaf.
 */
class class_test_set_linear_id: public TestDFS {
  void
  check_element () override
  {
    t8_debugf ("maxlvl: %i\n", maxlvl);
    t8_debugf ("level: %d\n", scheme->element_get_level (eclass, element));
    if (scheme->element_get_level (eclass, element) == maxlvl) {
      scheme->element_copy (eclass, element, test_element);
      scheme->element_set_linear_id (eclass, test_element, maxlvl, id_counter);
      id_counter++;
      EXPECT_ELEM_EQ (scheme, eclass, element, test_element);
    }
    else {
      return;
    }
  }

 protected:
  int
  set_maxlvl ()
  {
#if T8_ENABLE_LESS_TESTS
    const int maxlevel = 3;
#else
    const int maxlevel = 5;
#endif
    t8_debugf ("start maxlvl: %i\n", maxlevel);
    return maxlevel;
  }

  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    scheme->element_new (eclass, 1, &test_element);
    scheme->element_copy (eclass, element, test_element);
    maxlvl = set_maxlvl ();
  }
  void
  TearDown () override
  {
    /* Destroy element */
    scheme->element_destroy (eclass, 1, &test_element);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t *test_element;
  int id_counter = 0;
  int maxlvl;
};

TEST_P (class_test_set_linear_id, test_linear_id_dfs)
{
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_set_linear_id, AllSchemes, print_all_schemes);
