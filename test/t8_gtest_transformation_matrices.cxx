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

/**
 * In this file we collect test for projection matrices
 */

#include <gtest/gtest.h>
#include <t8_transformation_matrix.hxx>

#define epsilon 1e-9
typedef double      t8_test_point[3];

/**
 * Test the inverse camera transformation
 * The camera is centered at (0,0,0), facing points along the axis with dist 1. 
 * All points should be mapped to (0,0,-1) by the transformation.
 * 
 */
TEST (t8_gtest_transformation, cam_transformation)
{
  const t8_test_point cam = { 0.0, 0.0, 0.0 };
  /* Reference points on the axis */
  const t8_test_point ref_point[3] = { {1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0}
  };
  /* Upwards direction of the cam */
  const t8_test_point cam_up[3] = { {0.0, 0.0, 1.0},
  {0.0, 0.0, 1.0},
  {1.0, 0.0, 0.0}
  };
  /* All reference points should be mapped to this point */
  const t8_test_point result = { 0.0, 0.0, -1.0 };

  for (int i = 0; i < 3; i++) {
    t8_test_point       controll;

    double              cam_transform[4][4] = { 0.0 };
    /* Compute the transformation */
    inverse_camera_transformation (cam, ref_point[i], cam_up[i],
                                   cam_transform);
    mat4d_vec_multi (cam_transform, ref_point[i], controll);

    EXPECT_NEAR (result[0], controll[0], epsilon);
    EXPECT_NEAR (result[1], controll[1], epsilon);
    EXPECT_NEAR (result[2], controll[2], epsilon);
  }
}

/**
 * Test the perspective projection by checking if all the corners of the 
 * view-volume are mapped onto the corners of the translated cube [0,1]^3.
 */
TEST (t8_gtest_transformation, perspective_projection)
{
  /* Define what is going to be mapped. */
  /*The image width and the 0-point form a triangle with cathetes of equal length. */
  const double        near = 1.0;
  const double        far = 2.0 * near;
  const double        height = 1.0;
  const double        width = 1.0;
  double              perspective[4][4] = { 0.0 };

  /*The corners of the volume of view. */
  const t8_test_point test[8] = { {width, height, near},
  {-width, height, near},
  {width, -height, near},
  {-width, -height, near},
  {2 * width, 2 * height, far},
  {-2 * width, 2 * height, far},
  {2 * width, -2 * height, far},
  {-2 * width, -2 * height, far}
  };

  /* The coordinates of the projected test-points */
  const t8_test_point result[8] = { {-1.0, -1.0, 4.0},
  {1.0, -1.0, 4.0},
  {-1.0, 1.0, 4.0},
  {1.0, 1.0, 4.0},
  {-1.0, -1.0, 3.0},
  {1.0, -1.0, 3.0},
  {-1.0, 1.0, 3.0},
  {1.0, 1.0, 3.0}
  };

  perspective_projection (width, height, near, far, perspective);

  for (int i = 0; i < 8; i++) {
    t8_test_point       controll;
    mat4d_vec_multi (perspective, test[i], controll);
    EXPECT_NEAR (controll[0], result[i][0], epsilon);
    EXPECT_NEAR (controll[1], result[i][1], epsilon);
    EXPECT_NEAR (controll[2], result[i][2], epsilon);
  }
}
