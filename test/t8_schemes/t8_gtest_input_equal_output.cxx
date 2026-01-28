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
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"
#include <t8_data/t8_containers.h>

/* The mainly tested function in this test is element_get_children_at_face. The function allows that the input elem is equal children[0] (output).
Therefore the test checks if this is true and the input is not overridden in the loop or is overridden in the last iteration.
This is tested by comparing the output of two cases. In the first case the condition is not true in the second case elem is equal children[0]. */
struct class_test_equal: public TestDFS
{
 private:
  void
  check_element () override
  {
    for (int iface = 0; iface < scheme->element_get_num_faces (eclass, element); iface++) {

      const int num_children = scheme->element_get_num_face_children (eclass, element, iface);

      elems1 = T8_TESTSUITE_ALLOC (t8_element_t *, num_children);
      elems2 = T8_TESTSUITE_ALLOC (t8_element_t *, num_children);

      scheme->element_new (eclass, num_children, elems1);
      scheme->element_new (eclass, num_children, elems2);

      /* Copy element */
      t8_element_t *first_of_elems1 = elems1[0];
      scheme->element_copy (eclass, element, first_of_elems1);

      scheme->element_get_children_at_face (eclass, element, iface, elems2, num_children, NULL);

      scheme->element_get_children_at_face (eclass, first_of_elems1, iface, elems1, num_children, NULL);

      for (int ichild = 0; ichild < num_children; ichild++) {
        t8_element_t *child1 = elems1[ichild];
        t8_element_t *child2 = elems2[ichild];

        EXPECT_ELEM_EQ (scheme, eclass, child1, child2);
      }

      /* Destroy element */
      scheme->element_destroy (eclass, num_children, elems1);
      scheme->element_destroy (eclass, num_children, elems2);

      T8_TESTSUITE_FREE (elems1);
      T8_TESTSUITE_FREE (elems2);
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
  }
  void
  TearDown () override
  {
    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t **elems1;
  t8_element_t **elems2;
};

TEST_P (class_test_equal, test_equal_dfs)
{
#if T8_TEST_LEVEL_INT >= 1
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_equal, AllSchemes, print_all_schemes);
