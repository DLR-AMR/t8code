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

/** This test is used to test the corner_face and face_corner scheme functions.
 * The first test gets the corners of all faces of an element and checks if the reverse function can get the correct face of the corner.
 * The second test does the same but in reverse order.
 */
class class_face_corner_test: public TestDFS {

  void
  check_element () override
  {
    const int num_faces = scheme->element_get_num_faces (eclass, element);

    const int max_num_faces = scheme->element_get_max_num_faces (eclass, element);
    EXPECT_LE (num_faces, max_num_faces);

    for (int iface = 0; iface < num_faces; iface++) {
      const t8_element_shape_t face_shape = scheme->element_get_face_shape (eclass, element, iface);
      const int num_vertices = t8_element_shape_num_vertices (face_shape);
      for (int ivertex = 0; ivertex < num_vertices; ivertex++) {
        const int corner = scheme->element_get_face_corner (eclass, element, iface, ivertex);
        /* Only valid for hypercubes*/
        const int num_corner_faces = t8_eclass_to_dimension[eclass];

        bool found_face = false;
        for (int icorner_face = 0; icorner_face < num_corner_faces; icorner_face++) {
          const int face = scheme->element_get_corner_face (eclass, element, corner, icorner_face);
          if (face == iface) {
            found_face = true;
          }
        }
        EXPECT_TRUE (found_face);
      }
    }

    const int num_corners = scheme->element_get_num_corners (eclass, element);

    for (int icorner = 0; icorner < num_corners; icorner++) {
      /* Only valid for hypercubes*/
      const int num_faces = t8_eclass_to_dimension[eclass];
      for (int iface = 0; iface < num_faces; iface++) {
        const int face = scheme->element_get_corner_face (eclass, element, icorner, iface);

        const t8_element_shape_t face_shape = scheme->element_get_face_shape (eclass, element, face);
        const int num_face_corners = t8_element_shape_num_vertices (face_shape);

        bool found_corner = false;
        for (int iface_corner = 0; iface_corner < num_face_corners; iface_corner++) {
          const int corner = scheme->element_get_face_corner (eclass, element, face, iface_corner);
          if (corner == icorner) {
            found_corner = true;
          }
        }
        EXPECT_TRUE (found_corner);
      }
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
  }
  void
  TearDown () override
  {
    /* Destroy DFS test */
    dfs_test_teardown ();
  }
};

TEST_P (class_face_corner_test, test_equal_dfs)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_face_corner_test,
                          testing::Combine (AllSchemeCollections,
                                            ::testing::Range (T8_ECLASS_VERTEX, (t8_eclass_t) (T8_ECLASS_HEX + 1))),
                          print_all_schemes);
