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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include "t8_gtest_dfs_base.hxx"
#include <test/t8_gtest_macros.hxx>

class class_child_parent_face: public TestDFS {
  virtual void
  check_element ()
  {
    const int num_faces = ts->t8_element_num_faces (element);
    for (int iface = 0; iface < num_faces; iface++) {
      /* Iterate over all faces and determine the facechildren*/
      const int num_face_children = ts->t8_element_num_face_children (element, iface);
      t8_element_t **children;
      children = T8_ALLOC (t8_element_t *, num_face_children);
      ts->t8_element_new (num_face_children, children);

      ts->t8_element_children_at_face (element, iface, children, num_face_children, NULL);

      for (int ifacechild = 0; ifacechild < num_face_children; ifacechild++) {
        /* Iterate over those children and determine the childface corresponding to the parentface */
        int childface = ts->t8_element_face_child_face (element, iface, ifacechild);
        ASSERT_NE (childface, -1);
        /* Determine the parentface corresponding to the childface */
        int parentface = ts->t8_element_face_parent_face (children[ifacechild], childface);
        /* Check, that this is equal to the face that we started with */
        EXPECT_EQ (iface, parentface);
      }
      ts->t8_element_destroy (num_face_children, children);
      T8_FREE (children);
    }
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

TEST_P (class_child_parent_face, t8_recursive_dfs_child_parent_face)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 4;
#else
  const int maxlvl = 6;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_child_parent_face, class_child_parent_face, AllEclasses, print_eclass);
