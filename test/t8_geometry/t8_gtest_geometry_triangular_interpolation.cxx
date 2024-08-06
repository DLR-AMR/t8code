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
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry_helpers.h>

TEST_F (triangular_interpolation, corners_map_to_corners){
  const double eps = 1e-14;

  const double coeff_A[3] = {0, 0};
  const double coeff_B[3] = {0, 1};
  const double coeff_C[3] = {1, 1};
  const double coeff_D[3] = {1, 1};

  const double corner_values[] =
  {
    1, 2, // corner A
    3, 4, // corner B
    5, 6  // corner C
  };

  const double *cornerA = corner_values;
  const double *cornerB = corner_values + 2;
  const double *cornerC = corner_values + 4;

  double result[2];

  // Check corner A
  t8_geom_triangular_interpolation (coeff_A, corner_values, 2, 2, result);
  EXPECT_NEAR (cornerA[0], result[0], eps);
  EXPECT_NEAR (cornerA[1], result[1], eps);

  // Check corner B
  t8_geom_triangular_interpolation (coeff_B, corner_values, 2, 2, result);
  EXPECT_NEAR (cornerB[0], result[0], eps);
  EXPECT_NEAR (cornerB[1], result[1], eps);

  // Check corner C
  t8_geom_triangular_interpolation (coeff_C, corner_values, 2, 2, result);
  EXPECT_NEAR (cornerC[0], result[0], eps);
  EXPECT_NEAR (cornerC[1], result[1], eps);
}