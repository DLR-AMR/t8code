/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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
#include <t8_schemes/t8_consecutive/t8_consecutive_cxx.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include "t8_gtest_dfs_base.hxx"
#include <test/t8_gtest_macros.hxx>

class class_linear_id: public TestDFS {
  virtual void
  check_element ()
  {
    t8_debugf ("check element: \n");
    ts->t8_element_debug_print (element);
    linear_id = ts->t8_element_get_linear_id (element, ts->t8_element_level (element));
    t8_debugf ("linearid %li \n", linear_id);
    ts->t8_element_set_linear_id (test_element, ts->t8_element_level (element), linear_id);
    EXPECT_ELEM_EQ (ts, element, test_element);
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    ts->t8_element_new (1, &test_element);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    ts->t8_element_destroy (1, &test_element);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_gloidx_t linear_id;
  t8_element_t *test_element;
};

TEST_P (class_linear_id, t8_gtest_linear_id)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 2;
#else
  const int maxlvl = 2;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_linear_id, class_linear_id, AllEclasses);
