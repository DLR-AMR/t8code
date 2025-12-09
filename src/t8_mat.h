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

/** \file t8_mat.h
 * Collection of useful 3x3 matrices, matrix-matrix and matrix-vector operations.
 */

#ifndef T8_MAT_H
#define T8_MAT_H

#include <math.h>

/** Initialize given 3x3 matrix as rotation matrix around the x-axis with given angle.
 * \param [in,out]      mat     3x3-matrix.
 * \param [in]          angle   Rotation angle in radians.
 */
static inline void
t8_mat_init_xrot (double mat[3][3], const double angle)
{
  /* first row */
  mat[0][0] = 1.0;
  mat[0][1] = 0.0;
  mat[0][2] = 0.0;

  /* second row */
  mat[1][0] = 0.0;
  mat[1][1] = cos (angle);
  mat[1][2] = -sin (angle);

  /* third row */
  mat[2][0] = 0.0;
  mat[2][1] = sin (angle);
  mat[2][2] = cos (angle);
}

/** Initialize given 3x3 matrix as rotation matrix around the y-axis with given angle.
 * \param [in,out]      mat     3x3-matrix.
 * \param [in]          angle   Rotation angle in radians.
 */
static inline void
t8_mat_init_yrot (double mat[3][3], const double angle)
{
  /* first row */
  mat[0][0] = cos (angle);
  mat[0][1] = 0.0;
  mat[0][2] = sin (angle);

  /* second row */
  mat[1][0] = 0.0;
  mat[1][1] = 1.0;
  mat[1][2] = 0.0;

  /* third row */
  mat[2][0] = -sin (angle);
  mat[2][1] = 0.0;
  mat[2][2] = cos (angle);
}

/** Initialize given 3x3 matrix as rotation matrix around the z-axis with given angle.
 * \param [in,out]      mat     3x3-matrix.
 * \param [in]          angle   Rotation angle in radians.
 */
static inline void
t8_mat_init_zrot (double mat[3][3], const double angle)
{
  /* first row */
  mat[0][0] = cos (angle);
  mat[0][1] = -sin (angle);
  mat[0][2] = 0.0;

  /* second row */
  mat[1][0] = sin (angle);
  mat[1][1] = cos (angle);
  mat[1][2] = 0.0;

  /* third row */
  mat[2][0] = 0.0;
  mat[2][1] = 0.0;
  mat[2][2] = 1.0;
}

/** Apply matrix-matrix multiplication: b = M*a.
 * \param [in]        mat   3x3-matrix.
 * \param [in]        a     3-vector.
 * \param [in,out]    b     3-vector.
 */
static inline void
t8_mat_mult_vec (const double mat[3][3], const double a[3], double b[3])
{
  b[0] = mat[0][0] * a[0] + mat[0][1] * a[1] + mat[0][2] * a[2];
  b[1] = mat[1][0] * a[0] + mat[1][1] * a[1] + mat[1][2] * a[2];
  b[2] = mat[2][0] * a[0] + mat[2][1] * a[1] + mat[2][2] * a[2];
}

/** Apply matrix-matrix multiplication: C = A*B.
 * \param [in]        A     3x3-matrix.
 * \param [in]        B     3x3-matrix.
 * \param [in,out]    C     3x3-matrix.
 */
static inline void
t8_mat_mult_mat (const double A[3][3], const double B[3][3], double C[3][3])
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        C[i][j] = C[i][j] + A[i][k] * B[k][j];
      }
    }
  }
}

#endif /* !T8_MAT_H */
