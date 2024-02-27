/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_gtest_geometry_linear.cxx
 * This file contains tests for the linear axis aligned geometry of t8code.
 */

#include <test/t8_gtest_macros.hxx>
#include <test/t8_geometry/t8_gtest_geometry_macros.hxx>
#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_element.h>
#include <vector>

t8_cmesh_t
t8_geometry_testing_tree_from_class (const t8_eclass_t eclass)
{
  t8_cmesh_t cmesh;
  t8_geometry_c *geometry = new t8_geometry_linear_axis_aligned (t8_eclass_to_dimension[eclass]);
  t8_cmesh_init (&cmesh);
  std::vector<double> vertices (2 * T8_ECLASS_MAX_DIM, 0.0);
  std::fill (vertices.begin () + 3, vertices.end () - (T8_ECLASS_MAX_DIM - t8_eclass_to_dimension[eclass]), 1.0);
  t8_cmesh_set_tree_class (cmesh, 0, eclass);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices.data (), t8_eclass_num_vertices[eclass]);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  return cmesh;
}

/* Check whether the linear axis aligned geometry map is correct.
 * We create a cmesh of one tree and unit line, quad or hex
 * geometry. We then create random points in a reference tree and
 * check whether the evaluation is correct. */
TEST_P (geometry_test, cmesh_geometry_linear_axis_aligned)
{
  /* TODO: Add a test for the jacobian, as soon as its implemented. */

  t8_debugf ("Testing linear axis-aligned geometry evaluation.\n");

  /* Create random points in [0,1]^d and check if they are mapped correctly. */
  t8_geometry_linear_axis_aligned linear_axis_aligned_geom (t8_eclass_to_dimension[eclass]);
  t8_cmesh_t cmesh;

  double point[3];
  double point_mapped[3];
  const int seed = 0; /* RNG seed */
  const t8_geometry_c *cmesh_geom;

  cmesh = t8_geometry_testing_tree_from_class (eclass);

  /* Double check that the geometry is the linear geometry. */
  cmesh_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_TRUE (cmesh_geom != NULL) << "Could not get cmesh's geometry.";
  int has_same_name = strcmp (cmesh_geom->t8_geom_get_name (), linear_axis_aligned_geom.t8_geom_get_name ());
  ASSERT_EQ (has_same_name, 0) << "cmesh's geometry is not the linear geometry.";

  srand (seed);
  for (int ipoint = 0; ipoint < num_points; ++ipoint) {
    /* Compute random coordinates in [0,1].
     * These are seen as reference coordinates in the single
     * cmesh tree. Our geometry will map them into the physical
     * space. Since this space is also [0,1] and the cmesh only
     * has one tree, the mapped coordinates must be the same as the 
     * reference coordinates. */
    point[0] = (double) rand () / RAND_MAX;
    point[1] = (double) rand () / RAND_MAX;
    point[2] = (double) rand () / RAND_MAX;

    /* Evaluate the geometry */
    t8_geometry_evaluate (cmesh, 0, point, 1, point_mapped);
    /* Check that the first dim coordinates are the same */
    int idim;
    for (idim = 0; idim < t8_eclass_to_dimension[eclass]; ++idim) {
      ASSERT_NEAR (point[idim], point_mapped[idim], T8_PRECISION_SQRT_EPS)
        << "Linear axis aligned geometry computed wrong value.";
    }
    /* Check that the remaining entries are 0. */
    for (; idim < 3; ++idim) {
      ASSERT_EQ (point_mapped[idim], 0) << "Linear axis aligned geometry computed wrong value.";
    }
  }
  /* Destroy the cmesh */
  t8_cmesh_destroy (&cmesh);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry, geometry_test,
                          testing::Values (T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_HEX));
