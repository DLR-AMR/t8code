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

// typedef std::array<double, 3> t8_point_t;

template <std::size_t N>
using t8_point_t = std::array<double, N>;

#include <algorithm>
#include <numeric>
T8_EXTERN_C_BEGIN ();

/** Vector norm.
 * \param [in] vec  A 3D vector.
 * \return          The norm of \a vec.
 */
template <std::size_t N>
static inline double
t8_vec_norm (const t8_point_t<N> &vec)
{
  return std::sqrt (std::inner_product (vec.begin (), vec.end (), vec.begin (), 0.0));
}

/** Normalize a vector.
 * \param [in,out] vec  A 3D vector.
 */
template <std::size_t N>
static inline void
t8_vec_normalize (t8_point_t<N> vec)
{
  const double norm = t8_vec_norm (vec);
  std::transform (vec.begin (), vec.end (), vec.begin (), [norm] (double v) { return v / norm; });
}

/** Make a copy of a vector.
 * \param [in]  vec_in
 * \param [out] vec_out
 */
template <std::size_t N>
static inline void
t8_vec_copy (const t8_point_t<N> &vec_in, t8_point_t &vec_out)
{
  std::copy (vec_in.begin (), vec_in.end (), vec_out.begin ());
}

/** Euclidean distance of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The euclidean distance.
 *                     Equivalent to norm (X-Y).
 */
template <std::size_t N>
static inline double
t8_vec_dist (const t8_point_t<N> &vec_x, const t8_point_t<N> &vec_y)
{
  double dist = std::inner_product (vec_x.begin (), vec_x.end (), vec_y.begin (), 0.0, std::plus<double> (),
                                    [] (double x, double y) { return (x - y) * (x - y); });
  return std::sqrt (dist);
}

/** Compute X = alpha * X
 * \param [in,out] vec_x  A 3D vector. On output set to \a alpha * \a vec_x.
 * \param [in]     alpha  A factor.
 */
template <std::size_t N>
static inline void
t8_vec_ax (t8_point_t<N> &vec_x, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_x.begin (), [alpha] (double v) { return v * alpha; });
}

/** Compute Y = alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_z  On output set to \a alpha * \a vec_x.
 * \param [in]  alpha  A factor.
 */
template <std::size_t N>
static inline void
t8_vec_axy (const t8_point_t<N> vec_x, t8_point_t<N> vec_y, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), [alpha] (double v) { return v * alpha; });
}

/** Y = alpha * X + b
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_y  On input, a 3D vector.
 *                     On output set to \a alpha * \a vec_x + \a b.
 * \param [in]  alpha  A factor.
 * \param [in]  b      An offset.
 * \note It is possible that vec_x = vec_y on input to overwrite x
 */
template <std::size_t N>
static inline void
t8_vec_axb (const t8_point_t<N> vec_x, t8_point_t<N> vec_y, const double alpha, const double b)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), [alpha, b] (double v) { return alpha * v + b; });
}

/** Y = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in,out] vec_y On input, a 3D vector.
 *                      On output set \a to vec_y + \a alpha * \a vec_x
 * \param [in]  alpha  A factor.
 */
template <std::size_t N>
static inline void
t8_vec_axpy (const t8_point_t<N> vec_x, t8_point_t<N> vec_y, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), vec_y.begin (),
                  [alpha] (double x, double y) { return y + alpha * x; });
}

/** Z = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
 */
template <std::size_t N>
static inline void
t8_vec_axpyz (const t8_point_t<N> &vec_x, const t8_point_t<N> &vec_y, t8_point_t<N> &vec_z, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), vec_z.begin (),
                  [alpha] (double x, double y) { return y + alpha * x; });
}

/** Dot product of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The dot product \a vec_x * \a vec_y
 */
template <std::size_t N>
static inline double
t8_vec_dot (const t8_point_t<N> &vec_x, const t8_point_t<N> &vec_y)
{
  return std::inner_product (vec_x.begin (), vec_x.end (), vec_y.begin (), 0.0);
}
/** Cross product of X and Y
 * \param [in]  vec_x  A 2D vector.
 * \param [in]  vec_y  A 2D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
template <std::size_t N>
static inline void
t8_vec_cross_2D (const t8_point_t<N> vec_x, const t8_point_t<N> vec_y, double<N> &cross)
{
  cross = vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Cross product of X and Y
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
template <std::size_t N>
static inline void
t8_vec_cross_3D (const t8_point_t<N> vec_x, const t8_point_t<N> vec_y, t8_point_t<N> &cross)
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
template <std::size_t N>
static inline void
t8_vec_diff (const t8_point_t<N> vec_x, const t8_point_t<N> vec_y, double diff[3])
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), diff, std::minus<double> ());
}

/**
 * Check the equality of two vectors elementwise 
 * 
 * \param[in] vec_x 
 * \param[in] vec_y 
 * \param[in] tol 
 * \return true, if the vectors are equal up to \a tol 
 */
template <std::size_t N>
static inline bool
t8_vec_eq (const t8_point_t<N> &vec_x, const t8_point_t<N> &vec_y, const double tol)
{
  return std::equal (vec_x.begin (), vec_x.end (), vec_y.begin (),
                     [tol] (double x, double y) { return std::fabs (x - y) <= tol; });
}

/** Rescale a vector to a new length.
 * \param [in,out] vec  A 3D vector.
 * \param [in]  new_length  New length of the vector.
 */
template <std::size_t N>
static inline void
t8_vec_rescale (t8_point_t<N> &vec, const double new_length)
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
template <std::size_t N>
static inline void
t8_vec_tri_normal (const t8_point_t<N> p1, const t8_point_t<N> p2, const t8_point_t<N> p3, t8_point_t<N> normal)
{
  t8_point_t a;
  t8_point_t b;
  std::transform (p2.begin (), p2.end (), p1.begin (), a, std::minus<double> ());
  std::transform (p3.begin (), p3.end (), p1.begin (), b, std::minus<double> ());
  t8_vec_cross (a, b, normal);
}

/** Compute an orthogonal coordinate system from a given vector.
 * \param [in]   v1 3D vector.
 * \param [out]  v2 3D vector.
 * \param [out]  v3 3D vector.
 */
template <std::size_t N>
static inline void
t8_vec_orthogonal_tripod (const t8_point_t<N> &v1, t8_point_t<N> &v2, t8_point_t<N> &v3)
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
template <std::size_t N>
static inline void
t8_vec_swap (t8_point_t<N> &p1, t8_point_t<N> &p2)
{
  std::swap (p1, p2);
}

T8_EXTERN_C_END ();

#endif /* !T8_VEC_HXX */
