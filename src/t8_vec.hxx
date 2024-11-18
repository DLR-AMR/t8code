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

/** \file t8_vec.hxx
 * We define routines to handle 3-dimensional vectors.
 */

#ifndef T8_VEC_HXX
#define T8_VEC_HXX

#include <t8.h>

#include <algorithm>
#include <numeric>
T8_EXTERN_C_BEGIN ();

/** Vector norm.
 * \param [in] vec  A 3D vector.
 * \return          The norm of \a vec.
 */
static inline double
t8_vec_norm (const double vec[3])
{
  return std::sqrt (std::inner_product (vec, vec + 3, vec, 0.0));
}

/** Normalize a vector.
 * \param [in,out] vec  A 3D vector.
 */
static inline void
t8_vec_normalize (double vec[3])
{
  const double norm = t8_vec_norm (vec);
  std::transform (vec, vec + 3, vec, [norm] (double v) { return v / norm; });
}

/** Make a copy of a vector.
 * \param [in]  vec_in
 * \param [out] vec_out
 */
static inline void
t8_vec_copy (const double vec_in[3], double vec_out[3])
{
  std::copy (vec_in, vec_in + 3, vec_out);
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
  double dist = std::inner_product (vec_x, vec_x + 3, vec_y, 0.0, std::plus<double> (),
                                    [] (double x, double y) { return (x - y) * (x - y); });
  return std::sqrt (dist);
}

/** Compute X = alpha * X
 * \param [in,out] vec_x  A 3D vector. On output set to \a alpha * \a vec_x.
 * \param [in]     alpha  A factor.
 */
static inline void
t8_vec_ax (double vec_x[3], const double alpha)
{
  std::transform (vec_x, vec_x + 3, vec_x, [alpha] (double v) { return v * alpha; });
}

/** Compute Y = alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_z  On output set to \a alpha * \a vec_x.
 * \param [in]  alpha  A factor.
 */
static inline void
t8_vec_axy (const double vec_x[3], double vec_y[3], const double alpha)
{
  std::transform (vec_x, vec_x + 3, vec_y, [alpha] (double v) { return v * alpha; });
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
  std::transform (vec_x, vec_x + 3, vec_y, [alpha, b] (double v) { return alpha * v + b; });
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
  std::transform (vec_x, vec_x + 3, vec_y, vec_y, [alpha] (double x, double y) { return y + alpha * x; });
}

/** Z = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
 */
static inline void
t8_vec_axpyz (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha)
{
  std::transform (vec_x, vec_x + 3, vec_y, vec_z, [alpha] (double x, double y) { return y + alpha * x; });
}

/** Dot product of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The dot product \a vec_x * \a vec_y
 */
static inline double
t8_vec_dot (const double vec_x[3], const double vec_y[3])
{
  return std::inner_product (vec_x, vec_x + 3, vec_y, 0.0);
}

/** Cross product of X and Y
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
static inline void
t8_vec_cross (const double vec_x[3], const double vec_y[3], double cross[3])
{
  cross[0] = vec_x[1] * vec_y[2] - vec_x[2] * vec_y[1];
  cross[1] = vec_x[2] * vec_y[0] - vec_x[0] * vec_y[2];
  cross[2] = vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Compute the difference of two vectors.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] diff   On output, the difference of \a vec_x and \a vec_y.
 */
static inline void
t8_vec_diff (const double vec_x[3], const double vec_y[3], double diff[3])
{
  std::transform (vec_x, vec_x + 3, vec_y, diff, std::minus<double> ());
}

/**
 * Check the equality of two vectors elementwise 
 * 
 * \param[in] vec_x 
 * \param[in] vec_y 
 * \param[in] tol 
 * \return true, if the vectors are equal up to \a tol 
 */
static inline bool
t8_vec_eq (const double vec_x[3], const double vec_y[3], const double tol)
{
  return std::equal (vec_x, vec_x + 3, vec_y, [tol] (double x, double y) { return std::fabs (x - y) <= tol; });
}

/** Rescale a vector to a new length.
 * \param [in,out] vec  A 3D vector.
 * \param [in]  new_length  New length of the vector.
 */
static inline void
t8_vec_rescale (double vec[3], const double new_length)
{
  t8_vec_normalize (vec);
  t8_vec_ax (vec, new_length);
}

/** Compute the normal of a triangle given by its three vertices.
 * \param [in]  p1  A 3D vector.
 * \param [in]  p2  A 3D vector.
 * \param [in]  p3  A 3D vector.
 * \param [out] Normal vector of the triangle. (Not necessarily of length 1!)
 */
static inline void
t8_vec_tri_normal (const double p1[3], const double p2[3], const double p3[3], double normal[3])
{
  double a[3];
  double b[3];
  std::transform (p2, p2 + 3, p1, a, std::minus<double> ());
  std::transform (p3, p3 + 3, p1, b, std::minus<double> ());
  t8_vec_cross (a, b, normal);
}

/** Compute an orthogonal coordinate system from a given vector.
 * \param [in]   v1 3D vector.
 * \param [out]  v2 3D vector.
 * \param [out]  v3 3D vector.
 */
static inline void
t8_vec_orthogonal_tripod (const double v1[3], double v2[3], double v3[3])
{
  v2[0] = v1[1];
  v2[1] = v1[2];
  v2[2] = -v1[0];

  t8_vec_axpy (v1, v2, -t8_vec_dot (v1, v2));
  t8_vec_cross (v1, v2, v3);

  t8_vec_normalize (v2);
  t8_vec_normalize (v3);
}

/** Swap the components of two vectors.
 * \param [in,out]  p1  A 3D vector.
 * \param [in,out]  p2  A 3D vector.
 */
static inline void
t8_vec_swap (double p1[3], double p2[3])
{
  double tmp;
  for (int i = 0; i < 3; i++) {
    tmp = p1[i];
    p1[i] = p2[i];
    p2[i] = tmp;
  }
}

T8_EXTERN_C_END ();

#endif /* !T8_VEC_HXX */
