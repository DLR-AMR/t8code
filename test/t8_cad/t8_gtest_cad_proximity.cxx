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
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_cad/t8_cad_shape_proximity.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Cut.hxx>

/* In this file we collect tests for t8code's CAD proximity module.
 * These tests are
 *  - element_and_point_inside: Builds a forest and a box shape and checks if the elements 
 *                              and points of the forest are clssified correctly.
 *  - element_skewed:           Builds a forest, which is not axis aligned and checks if CAD
 *                              throws an error.
 */

/**
 * Creates a axis aligned box with a box cutout. The cutout is deformed in
 * such a way, that it can also test the tolerances of the algorithms.
 * 
 * \return               The box shape.
 */
TopoDS_Shape
t8_gtest_cad_proximity_create_box_with_cutout ()
{
  gp_Pnt              pnt1 (0.4, 0.4, 0.4), pnt2 (1, 1, 1);
  TopoDS_Shape        argument = BRepPrimAPI_MakeBox (pnt1, pnt2);
  TopoDS_Shape        tool = BRepPrimAPI_MakeBox (0.75,
                                                  0.75 -
                                                  Precision::Confusion () * 2,
                                                  0.75 +
                                                  Precision::Confusion () *
                                                  2);
  return BRepAlgoAPI_Cut (argument, tool);
}

/**
 * Tests if a coordinate lies inside the box created in 
 * \a t8_gtest_cad_proximity_create_box_with_cutout.
 * Respects the given tolerances. Only defined in [0, 1]^3.
 * 
 * \param [in] coords    The coordinates to test. 
 * \return               True if point is inside the box.
 */
int
t8_gtest_cad_proximity_check_coordinate (double *coords)
{
  T8_ASSERT ((0 <= coords[0] && coords[0] <= 1) ||
             (0 <= coords[1] && coords[1] <= 1) ||
             (0 <= coords[2] && coords[2] <= 1));

  if (coords[0] <= 0.4 - Precision::Confusion () ||
      coords[1] <= 0.4 - Precision::Confusion () ||
      coords[2] <= 0.4 - Precision::Confusion ()) {
    return 0;
  }
  else if (coords[0] >= 0.75 ||
           coords[1] >= 0.75 - Precision::Confusion () * 2 ||
           coords[2] >= 0.75 + Precision::Confusion () * 2) {
    return 1;
  }
  else {
    return 0;
  }
}

/**
 * Creates a hypercube forest.
 * 
 * \param [in] level     The level of the forest.
 * \param [in] comm      The communicator to use.
 * \return               The forest.
 */
t8_forest_t
t8_gtest_cad_proximity_create_forest (int level, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_scheme_cxx_t    *scheme;
  scheme = t8_scheme_new_default_cxx ();
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  return t8_forest_new_uniform (cmesh, scheme, level, 0, comm);
}

/**
 * Tests the element and point inside functions of \a t8_cad_shape_proximity.
 */
TEST (t8_gtest_cad_proximity, element_and_point_inside)
{
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;
  t8_forest_t         forest;
  t8_element_t       *element;
  double              coords[3];
  TopoDS_Shape        shape;
  /* Generate the shape, cad object and forest. */
  shape = t8_gtest_cad_proximity_create_box_with_cutout ();
  t8_cad_shape_proximity cad (shape, 0);
#ifdef T8_ENABLE_LESS_TESTS
  forest = t8_gtest_cad_proximity_create_forest (4, comm);
#else
  forest = t8_gtest_cad_proximity_create_forest (2, comm);
#endif /* T8_ENABLE_LESS_TESTS */
  const t8_locidx_t   num_elements_in_tree =
    t8_forest_get_tree_num_elements (forest, 0);
  /* Iterate over each element of the forest and check,
   * if the element or point is inside. The box is designed in such a way,
   * that an element is inside the box, if its corner 7 is inside the box.
   * Therefore, we can check the point inside as well es the element inside
   * algorithm in the same loop. */
  for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
    element = t8_forest_get_element_in_tree (forest, 0, ielement);
    t8_forest_element_coordinate (forest, 0, element, 7, coords);
    const int           expect =
      t8_gtest_cad_proximity_check_coordinate (coords);
    const int           elem_inside =
      cad.t8_cad_is_element_inside_shape (forest, 0,
                                          element, 0, 1);
    const int           point_inside =
      (int) cad.t8_cad_is_point_inside_shape (coords, 0);
    EXPECT_EQ (elem_inside, expect);
    EXPECT_EQ (point_inside, expect);
  }
  t8_forest_unref (&forest);
}
