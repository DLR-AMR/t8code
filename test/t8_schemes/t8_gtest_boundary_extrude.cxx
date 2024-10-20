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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"
#include <test/t8_gtest_macros.hxx>

class class_test_boundary_extrude: public TestDFS {
  /* For elements that are on the face of the root element, check that creating the boundary element
   * and extruding it results in the original element
    */
  virtual void
  check_element ()
  {
    const int num_faces = ts->t8_element_num_faces (element);
    for (int iface = 0; iface < num_faces; iface++) {
      /* Iterate over all faces that are also root faces and determine the face element */
      if (ts->t8_element_is_root_boundary (element, iface)) {
        /* Get face scheme */
        const int tree_face = ts->t8_element_tree_face (element, iface);
        const t8_eclass_t face_eclass = (t8_eclass_t) t8_eclass_face_types[eclass][tree_face];
        const t8_scheme *face_ts = scheme->eclass_schemes[face_eclass];

        t8_element_t *boundary;
        face_ts->t8_element_new (1, &boundary);

        ts->t8_element_boundary_face (element, iface, boundary, face_ts);

        ts->t8_element_extrude_face (boundary, face_ts, check, tree_face);

        EXPECT_ELEM_EQ (ts, element, check);

        face_ts->t8_element_destroy (1, &boundary);
      }
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    /* Get element and initialize it */
    ts->t8_element_new (1, &check);
  }
  void
  TearDown () override
  {
    /* Destroy element */
    ts->t8_element_destroy (1, &check);

    /* Destroy DFS test */
    dfs_test_teardown ();
  }
  t8_element_t *check;
};

TEST_P (class_test_boundary_extrude, test_boundary_extrude_dfs)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 4;
#else
  const int maxlvl = 6;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_test_boundary_extrude, AllEclasses, print_eclass);
