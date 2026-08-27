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

/* In this file we collect tests for the routines in `t8_mat.h`. */

#include <t8.h>
#include <t8_helper_functions/t8_mat.h>
#include <gtest/gtest.h>

/* Test the `t8_mat_mult_vec` function. Here we also test the rotation matrices
 * at the same time. */
TEST (t8_gtest_mat, mat_mult_vec)
{
  const double start_vec[3] = { 0.0, 0.0, 1.0 };

  double xrotated_vec[3];
  double yrotated_vec[3];
  double zrotated_vec[3];

  /* Rotate -45 degrees around x-axis. */
  {

    /* Apply rotation. */
    double rot_mat[3][3];
    const double angle = -M_PI_4; /* -45 degrees */
    t8_mat_init_xrot (rot_mat, angle);
    t8_mat_mult_vec (rot_mat, start_vec, xrotated_vec);

    /* Check results. */
    const double test_vec[3] = { 0.0, 0.5 * M_SQRT2, 0.5 * M_SQRT2 };
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR (test_vec[i], xrotated_vec[i], T8_PRECISION_EPS);
    }
  }

  /* Rotate -45 degrees around z-axis. */
  {
    /* Apply rotation. */
    double rot_mat[3][3];
    const double angle = -M_PI_4; /* -45 degrees */
    t8_mat_init_zrot (rot_mat, angle);
    t8_mat_mult_vec (rot_mat, xrotated_vec, yrotated_vec);

    /* Check results. */
    const double test_vec[3] = { 0.5, 0.5, 0.5 * M_SQRT2 };
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR (test_vec[i], yrotated_vec[i], T8_PRECISION_EPS);
    }
  }

  /* Rotate 90 degrees around y-axis. */
  {
    /* Apply rotation. */
    double rot_mat[3][3];
    const double angle = M_PI_2; /* 90 degrees */
    t8_mat_init_yrot (rot_mat, angle);
    t8_mat_mult_vec (rot_mat, yrotated_vec, zrotated_vec);

    /* Check results. */
    const double test_vec[3] = { 0.5 * M_SQRT2, 0.5, -0.5 };
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR (test_vec[i], zrotated_vec[i], T8_PRECISION_EPS);
    }
  }
}

/* Test the `t8_mat_mult_mat` function. We replicate the three rotations from
 * above by cascading the matrix operators.*/
TEST (t8_gtest_mat, mat_mult_mat)
{
  const double start_vec[3] = { 0.0, 0.0, 1.0 };

  /* First rotate -45 degrees around x-axis. Then
   * rotate -45 degrees around z-axis.
   * Finally, rotate 90 degrees around the y-axis. */

  double xrot_mat[3][3];
  double yrot_mat[3][3];
  double zrot_mat[3][3];

  double xzrot_mat[3][3];
  double xzyrot_mat[3][3];

  t8_mat_init_xrot (xrot_mat, -M_PI_4);
  t8_mat_init_zrot (zrot_mat, -M_PI_4);
  t8_mat_init_yrot (yrot_mat, M_PI_2);

  t8_mat_mult_mat (zrot_mat, xrot_mat, xzrot_mat);
  t8_mat_mult_mat (yrot_mat, xzrot_mat, xzyrot_mat);

  {
    /* Apply rotation. */
    double rotated_vec[3];
    t8_mat_mult_vec (xzrot_mat, start_vec, rotated_vec);

    /* Check results. */
    const double test_vec[3] = { 0.5, 0.5, 0.5 * M_SQRT2 };
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR (test_vec[i], rotated_vec[i], T8_PRECISION_EPS);
    }
  }

  {
    /* Apply rotation. */
    double rotated_vec[3];
    t8_mat_mult_vec (xzyrot_mat, start_vec, rotated_vec);

    /* Check results. */
    const double test_vec[3] = { 0.5 * M_SQRT2, 0.5, -0.5 };
    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR (test_vec[i], rotated_vec[i], T8_PRECISION_EPS);
    }
  }
}
