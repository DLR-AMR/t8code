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

/** \file t8_vec.h
 * We define routines to handle 3-dimensional vectors.
 */

#ifndef T8_VEC_H
#define T8_VEC_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Vector norm.
 * \param [in] vec  A 3D vector.
 * \return          The norm of \a vec.
 */
static inline double
t8_vec_norm (const double vec[3])
{
  double norm = 0;

  for (int i = 0; i < 3; i++) {
    norm += vec[i] * vec[i];
  }
  return sqrt (norm);
}

/** Normalize a vector.
 * \param [in,out] vec  A 3D vector.
 */
static inline void
t8_vec_normalize (double vec[3])
{
  const double inv_norm = 1.0 / t8_vec_norm (vec);

  for (int i = 0; i < 3; i++) {
    vec[i] *= inv_norm;
  }
}

/** Euclidean distance of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The euclidean distance.
 *                     Equivalent to norm (X-Y).
 */
static inline double
t8_vec_dist (const double vec_x[3], const double vec_y[3])
{
  double dist = 0;

  for (int i = 0; i < 3; i++) {
    dist += SC_SQR (vec_x[i] - vec_y[i]);
  }
  return sqrt (dist);
}

/** Compute X = alpha * X
 * \param [in,out] vec_x  A 3D vector. On output set to \a alpha * \a vec_x.
 * \param [in]     alpha  A factor.
 */
static inline void
t8_vec_ax (double vec_x[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_x[i] *= alpha;
  }
}

/** Compute Y = alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_z  On output set to \a alpha * \a vec_x.
 * \param [in]  alpha  A factor.
 */
static inline void
t8_vec_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] = vec_x[i] * alpha;
  }
}

/** Y = alpha * X + b
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_y  On input, a 3D vector.
 *                     On output set to \a alpha * \a vec_x + \a b.
 * \param [in]  alpha  A factor.
 * \param [in]  b      An offset.
 * \note It is possible that vec_x = vec_y on input to overwrite x
 */
static inline void
t8_vec_axb (const double vec_x[3], double vec_y[3], const double alpha, const double b)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] = alpha * vec_x[i] + b;
  }
}

/** Y = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in,out] vec_y On input, a 3D vector.
 *                      On output set \a to vec_y + \a alpha * \a vec_x
 * \param [in]  alpha  A factor.
 */
static inline void
t8_vec_axpy (const double vec_x[3], double vec_y[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_y[i] += alpha * vec_x[i];
  }
}

/** Z = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
 */
static inline void
t8_vec_axpyz (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  for (int i = 0; i < 3; i++) {
    vec_z[i] = vec_y[i] + alpha * vec_x[i];
  }
}

/** Dot product of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The dot product \a vec_x * \a vec_y
 */
static inline double
t8_vec_dot (const double vec_x[3], const double vec_y[3])
{
  double dot = 0;

  for (int i = 0; i < 3; i++) {
    dot += vec_x[i] * vec_y[i];
  }
  return dot;
}

/** Cross product of X and Y
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
static inline void
t8_vec_cross (const double vec_x[3], const double vec_y[3], double cross[3])
{
  for (int i = 0; i < 3; i++) {
    cross[i] = vec_x[(i + 1) % 3] * vec_y[(i + 2) % 3] - vec_x[(i + 2) % 3] * vec_y[(i + 1) % 3];
  }
}

/** Compute the difference of two vectors.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] diff   On output, the difference of \a vec_x and \a vec_y.
 */
static inline void
t8_vec_diff (const double vec_x[3], const double vec_y[3], double diff[3])
{
  for (int i = 0; i < 3; i++) {
    diff[i] = vec_x[i] - vec_y[i];
  }
}

/**
 * Check the equality of two vectors elementwise 
 * 
 * \param[in] vec_x 
 * \param[in] vec_y 
 * \param[in] eps 
 * \return true, if the vectors are equal up to \a eps 
 */
static inline int
t8_vec_eq (const double vec_x[3], const double vec_y[3], const double eps)
{
  return fabs (vec_x[0] - vec_y[0]) < eps && fabs (vec_x[1] - vec_y[1]) < eps && fabs (vec_x[2] - vec_y[2]) < eps;
}

T8_EXTERN_C_END ();

#endif /* !T8_VEC_H! */
