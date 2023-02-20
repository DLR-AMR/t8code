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

#include<t8_transformation_matrix.hxx>
#include<t8_vec.h>

void
fill_mat4d (double mat[4][4], const double data[16])
{
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      mat[i][j] = data[4 * i + j];
    }
  }
}

void
mat4d_vec_multi (const double mat[4][4], const double vec[3],
                 double out_vec[3])
{
  out_vec[0] =
    mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3];
  out_vec[1] =
    mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3];
  out_vec[2] =
    mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + mat[2][3];
  const double        w =
    mat[3][0] * vec[0] + mat[3][1] * vec[1] + mat[3][2] * vec[2] + mat[3][3];

  /* From homogenous to cartesian coordinates */
  if (w != 1.0) {
    T8_ASSERT (w != 0.0);
    out_vec[0] /= w;
    out_vec[1] /= w;
    out_vec[2] /= w;
  }
}

void
inverse_camera_transformation (const double cam[3], const double ref_point[3],
                               const double up[3],
                               double transformation[4][4])
{
  double              cam_x_axis[3] = { 0.0 };
  double              cam_y_axis[3] = { 0.0 };
  double              cam_z_axis[3] = { cam[0], cam[1], cam[2] };
  /* Compute the vector from the view_point to the camera. */
  t8_vec_axpy (ref_point, cam_z_axis, -1.0);
  /* Compute the norm of the previously computed vector. */
  const double        cam_z_norm = t8_vec_norm (cam_z_axis);
  /* normalize the z-axis of the cam-coord-system. */
  t8_vec_ax (cam_z_axis, cam_z_norm);
  /* The x-axis of the camera system is the dot product of the normalized up-axis and
   * and the camera-z-axis. */
  const double        up_norm = t8_vec_norm (up);
  const double        up_normalized[3] =
    { up[0] / up_norm, up[1] / up_norm, up[2] / up_norm };
  t8_vec_cross (up_normalized, cam_z_axis, cam_x_axis);

  /* The y-axis of the camera. */
  t8_vec_cross (cam_z_axis, cam_x_axis, cam_y_axis);

  /*Set-up the transformation-matrix */
  const double        trans[3][3] =
    { {cam_x_axis[0], cam_x_axis[1], cam_x_axis[2]},
  {cam_y_axis[0], cam_y_axis[1], cam_y_axis[2]},
  {cam_z_axis[0], cam_z_axis[1], cam_z_axis[2]}
  };
  double              translation[3] = { 0.0 };
  t8_mat_vec (trans, cam, -1.0, translation);

  const double        data[16] =
    { cam_x_axis[0], cam_x_axis[1], cam_x_axis[2], translation[0],
    cam_y_axis[0], cam_y_axis[1], cam_y_axis[2], translation[1],
    cam_z_axis[0], cam_z_axis[1], cam_z_axis[2], translation[2],
    0.0, 0.0, 0.0, 1.0
  };
  fill_mat4d (transformation, data);
  return;
}

void
perspective_projection (const double width,
                        const double height, const double near,
                        const double far, double projection[4][4])
{
  /* near and width form a triangle with an angle of 90deg between near and width.
   * The field of view is the angle between near and the hypothenuse. */
  const double        hypothenuse = sqrt (near * near + width * width);
  const double        fov = asin (width / hypothenuse);

  const double        scale = cos (fov) / sin (fov);
  const double        aspect = width / height;
  const double        data[16] = { scale / aspect, 0.0, 0.0, 0.0,
    0.0, scale, 0.0, 0.0,
    0.0, 0.0, (-far) / (far - near), (-far * near) / (far - near),
    0.0, 0.0, -1.0, 0.0
  };
  fill_mat4d (projection, data);
}
