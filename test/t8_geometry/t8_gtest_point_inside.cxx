/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
#include <sc_functions.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_element_cxx.hxx>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>

/* In this test we define a triangle in the x-y plane
 * and a point that lies in a triangle that is parallel
 * to this triangle on the z-axis.
 * The point must be correctly identified as lying outside
 * of the triangle.
 */
TEST (t8_point_inside, test_point_inside_specific_triangle)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_element_t       *element;
  double              vertices[9] = {
    0., 0., 0.,
    1., 0., 0.,
    1., 1., 0.
  };
  double              test_point[3] = {
    0.3, 0.3, 1
  };
  int                 point_is_inside;
  const double        tolerance = 1e-12;        /* Numerical tolerance that we allow for the point inside check */
  t8_geometry_c      *linear_geom = new t8_geometry_linear (2);

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  /* We use standard linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);

  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0,
                           sc_MPI_COMM_WORLD);

  if (t8_forest_get_local_num_elements (forest) <= 0) {
    /* Skip empty forests (can occur when executed in parallel) */
    t8_forest_unref (&forest);
    GTEST_SKIP ();
  }

  element = t8_forest_get_element (forest, 0, NULL);

  point_is_inside =
    t8_forest_element_point_inside (forest, 0, element, test_point,
                                    tolerance);
  ASSERT_FALSE (point_is_inside) <<
    "The point is wrongly detected as inside the triangle.";
  t8_forest_unref (&forest);
}

/* In this test we define a quad in the x-y plane
 * and a point that lies in a quad that is parallel
 * to the first quad on the z-axis.
 * The point must be correctly identified as lying outside
 * of the quad.
 */
TEST (t8_point_inside, test_point_inside_specific_quad)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_element_t       *element;
  double              vertices[12] = {
    0., 0., 0.,
    1., 0., 0.,
    0., 1., 0.,
    1., 1., 0.
  };
  double              test_point[3] = {
    0.3, 0.3, 1
  };
  int                 point_is_inside;
  const double        tolerance = 1e-12;        /* Numerical tolerance that we allow for the point inside check */
  t8_geometry_c      *linear_geom = new t8_geometry_linear (2);

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  /* We use standard linear geometry */
  t8_cmesh_register_geometry (cmesh, linear_geom);

  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0,
                           sc_MPI_COMM_WORLD);

  if (t8_forest_get_local_num_elements (forest) <= 0) {
    /* Skip empty forests (can occur when executed in parallel) */
    t8_forest_unref (&forest);
    GTEST_SKIP ();
  }

  element = t8_forest_get_element (forest, 0, NULL);

  point_is_inside =
    t8_forest_element_point_inside (forest, 0,
                                    element, test_point, tolerance);

  ASSERT_FALSE (point_is_inside) <<
    "The point is wrongly detected as inside the triangle.";

  t8_forest_unref (&forest);
}
