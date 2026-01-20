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
#include <t8_eclass.h>
#include <t8_gtest_schemes.hxx>
#include <t8_gtest_custom_assertion.hxx>
#include <t8_gtest_macros.hxx>
#include "t8_gtest_bfs_base.hxx"

/** In this test we iterate through all elements. 
 * On every level we check if the element is equal to the element we get when setting it from the linear id.
 * The id_counter is then increased to match the id of the next leaf. After we have reached the last element on a level, 
 * we increase the level and reset the id_counter to 0.
 */
class class_test_set_linear_id: public TestBFS {
  void
  check_element () override
  {
    if (computed_level < current_level) {
      computed_level = current_level;
      id_counter = 0;
    }
    scheme->element_set_linear_id (eclass, test_element, current_level, id_counter);
    id_counter++;
    EXPECT_ELEM_EQ (scheme, eclass, element, test_element);
  }

 protected:
  void
  SetUp () override
  {
    bfs_test_setup ();
    /* Get element and initialize it */
    scheme->element_new (eclass, 1, &test_element);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    scheme->element_destroy (eclass, 1, &test_element);

    /* Destroy BFS test */
    bfs_test_teardown ();
  }
#if T8_TEST_LEVEL_INT >= 1
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif
  int computed_level = 0;
  t8_linearidx_t id_counter = 0;
  t8_element_t *test_element;
};

TEST_P (class_test_set_linear_id, test_linear_id_bfs)
{
  check_recursive_bfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_set_linear_id, AllSchemes, print_all_schemes);
