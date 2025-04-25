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
class class_reference_coords_test: public TestDFS {

  void
  check_element () override
  {
    for (int ivertex = 0; ivertex < scheme->element_get_num_corners (eclass, element); ivertex++) {
      double vertex_coords[3];
      double reference_coords[3];
      scheme->element_debug_print (eclass, element);
      scheme->element_get_vertex_reference_coords (eclass, element, ivertex, vertex_coords);
      t8_debugf ("vertex_coords: %lf, %lf, %lf \n", vertex_coords[0], vertex_coords[1], vertex_coords[2]);
      t8_element_shape_t shape = scheme->element_get_shape (eclass, element);
      scheme->element_get_reference_coords (eclass, element, t8_element_corner_ref_coords[shape][ivertex], 1,
                                            reference_coords);
      t8_debugf ("reference_lut: %lf, %lf, %lf \n", t8_element_corner_ref_coords[shape][ivertex][0],
                 t8_element_corner_ref_coords[shape][ivertex][1], t8_element_corner_ref_coords[shape][ivertex][2]);
      t8_debugf ("reference_coords: %lf, %lf, %lf \n", reference_coords[0], reference_coords[1], reference_coords[2]);

      for (int idim = 0; idim < T8_ELEMENT_DIM[eclass]; idim++) {
        EXPECT_NEAR (vertex_coords[idim], reference_coords[idim], T8_PRECISION_SQRT_EPS);
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

TEST_P (class_reference_coords_test, test_equal_dfs)
{
#if T8CODE_TEST_LEVEL >= 1
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif
  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_test_all_imps, class_reference_coords_test, AllSchemes, print_all_schemes);
