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

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry_helpers.h>

/* In this file we test the functionality of t8_geom_triangular_interpolation.
 * The function maps reference coordinates inside the reference triangle or tetrahedron
 * to the coordinate space of a triangle or tetrahedron given by its corner coordinates.
 * We currently check that the reference coordinates of the corners are mapped correctly
 * to the corners of the interpolated triangle/tet. */

/* TODO:
 *       - Check corner_value_dim > 1
 *       - Check different coefficient coordinates (for example midpoint)
 *       - Check different corner values
 *       - The 2d and 3d example contain a lot of nearly duplicate code.
 *         Can we simplify?
 */

/* Check that the corner interpolation coordinates are
 * mapped to the corners of a triangle. */
TEST (triangular_interpolation, corners_map_to_corners_2d)
{
  const double coeff_A[2] = { 0, 0 };
  const double coeff_B[2] = { 0, 1 };
  const double coeff_C[2] = { 1, 1 };

  const double corner_values_2d[] = {
    1, 2,  // corner A
    3, 4,  // corner B
    5, 6   // corner C
  };

  const double *cornerA = corner_values_2d;
  const double *cornerB = corner_values_2d + 2;
  const double *cornerC = corner_values_2d + 4;

  const double dimension = 2;

  double result[2];

  // Check corner A
  // 2D
  t8_geom_triangular_interpolation (coeff_A, corner_values_2d, 2, dimension, result);
  EXPECT_NEAR (cornerA[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerA[1], result[1], T8_PRECISION_EPS);

  // Check corner B
  // 2D
  t8_geom_triangular_interpolation (coeff_B, corner_values_2d, 2, dimension, result);
  EXPECT_NEAR (cornerB[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerB[1], result[1], T8_PRECISION_EPS);

  // Check corner C
  // 2D
  t8_geom_triangular_interpolation (coeff_C, corner_values_2d, 2, dimension, result);
  EXPECT_NEAR (cornerC[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerC[1], result[1], T8_PRECISION_EPS);
}

/* Check that the corner interpolation coordinates are
 * mapped to the corners of a tetrahedron. */
TEST (triangular_interpolation, corners_map_to_corners_3d)
{
  const double coeff_A[3] = { 0, 0, 0 };
  const double coeff_B[3] = { 0, 1, 0 };
  const double coeff_C[3] = { 1, 1, 0 };
  const double coeff_D[3] = { 1, 1, 1 };

  const double corner_values_3d[] = {
    1,  2,  3,   // corner A
    4,  5,  6,   // corner B
    7,  8,  9,   // corner C
    10, 11, 12,  // corner D
  };

  const double *cornerA = corner_values_3d;
  const double *cornerB = corner_values_3d + 3;
  const double *cornerC = corner_values_3d + 6;
  const double *cornerD = corner_values_3d + 9;

  const double dimension = 3;

  double result[3];

  // Check corner A
  // 3D
  t8_geom_triangular_interpolation (coeff_A, corner_values_3d, 3, dimension, result);
  EXPECT_NEAR (cornerA[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerA[1], result[1], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerA[2], result[2], T8_PRECISION_EPS);

  // Check corner B
  // 3D
  t8_geom_triangular_interpolation (coeff_B, corner_values_3d, 3, dimension, result);
  EXPECT_NEAR (cornerB[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerB[1], result[1], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerB[2], result[2], T8_PRECISION_EPS);

  // Check corner C
  // 3D
  t8_geom_triangular_interpolation (coeff_C, corner_values_3d, 3, dimension, result);
  EXPECT_NEAR (cornerC[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerC[1], result[1], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerC[2], result[2], T8_PRECISION_EPS);

  // Check corner D
  // 3D
  t8_geom_triangular_interpolation (coeff_D, corner_values_3d, 3, dimension, result);
  EXPECT_NEAR (cornerD[0], result[0], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerD[1], result[1], T8_PRECISION_EPS);
  EXPECT_NEAR (cornerD[2], result[2], T8_PRECISION_EPS);
}
