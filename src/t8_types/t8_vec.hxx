/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <t8_types/t8_type.hxx>
#include <t8_types/t8_operators.hxx>
#include <t8.h>
// template <std::size_t N>
// using t8_point_t = std::array<double, N>;
// typedef t8_point_t<3> t8_3D_point;

#include <algorithm>
#include <numeric>

template <std::size_t dim>
struct t8_vec_tag
{
};

template <std::size_t dim>
struct t8_point_tag
{
};

/**
 * Type alias for a vector.
 * \tparam dim Dimension of the vector.
 */
template <std::size_t dim>
using t8_vec = T8Type<std::array<double, dim>, t8_vec_tag<dim>, EqualityComparable, Swapable, RandomAccessible>;

/** Type alias for a 3D vector.
 * \note This is a convenience type alias for 3D vectors.
 */
using t8_3D_vec = t8_vec<3>;

/** 
 * Type alias for a point in N-dimensional space.
 * \tparam dim Dimension of the point.
 */
template <std::size_t dim>
using t8_point = T8Type<std::array<double, dim>, t8_point_tag<dim>, EqualityComparable, Swapable, RandomAccessible>;

/** Type alias for a 3D point.
 * \note This is a convenience type alias for 3D points.
 */
using t8_3D_point = t8_point<3>;

/** Vector norm.
 * \param [in] vec  An N-dimensional vector.
 * \return          The norm of \a vec.
 */
template <std::size_t dim>
static inline double
t8_norm (const t8_vec<dim> &vec)
{
  return std::sqrt (std::inner_product (vec.begin (), vec.end (), vec.begin (), 0.0));
}

/** Normalize a vector.
 * \param [in,out] vec  An N-dimensional vector.
 */
template <std::size_t dim>
constexpr void
t8_normalize (t8_vec<dim> &vec)
{
  const double norm = t8_norm (vec);
  std::transform (vec.begin (), vec.end (), vec.begin (), [norm] (double v) { return v / norm; });
}

/** Make a copy of a vector or point.
 * \param [in]  src  The source vector or point.
 * \param [out] dest The destination vector or point.
 */
template <typename T>
static inline void
t8_copy (const T &src, T &dest);

/**
 * Copy a vector.
 * \param [in]  src  The source vector.
 * \param [out] dest The destination vector.
 */
template <std::size_t dim>
constexpr void
t8_copy (const t8_vec<dim> &src, t8_vec<dim> &dest)
{
  std::copy (src.begin (), src.end (), dest.begin ());
}

/**
 * Copy a point.
 * \param [in]  src  The source point.
 * \param [out] dest The destination point.
 */
template <std::size_t dim>
constexpr void
t8_copy (const t8_point<dim> &src, t8_point<dim> &dest)
{
  std::copy (src.begin (), src.end (), dest.begin ());
}

/** Euclidean distance of X and Y.
 * \param [in]  point_x  An N-dimensional point.
 * \param [in]  point_y  An N-dimensional point.
 * \return             The euclidean distance.
 *                     Equivalent to norm (X-Y).
 */
template <std::size_t dim>
static inline double
t8_dist (const t8_point<dim> &point_x, const t8_point<dim> &point_y)
{
  double dist = std::inner_product (point_x.begin (), point_x.end (), point_y.begin (), 0.0, std::plus<double> (),
                                    [] (double x, double y) { return (x - y) * (x - y); });
  return std::sqrt (dist);
}

/** Compute X = alpha * X
 * \param [in,out] vec_x  An N-dimensional vector. On output set to \a alpha * \a vec_x.
 * \param [in]     alpha  A factor.
 */
template <std::size_t dim>
constexpr void
t8_ax (t8_vec<dim> &vec_x, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_x.begin (), [alpha] (double v) { return v * alpha; });
}

/** Compute Y = alpha * X
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [out] vec_y  On output set to \a alpha * \a vec_x.
 * \param [in]  alpha  A factor.
 */
template <std::size_t dim>
constexpr void
t8_axy (const t8_vec<dim> &vec_x, t8_vec<dim> &vec_y, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), [alpha] (double v) { return v * alpha; });
}

/** Y = alpha * X + b
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [out] vec_y  On input, An N-dimensional vector.
 *                     On output set to \a alpha * \a vec_x + \a b.
 * \param [in]  alpha  A factor.
 * \param [in]  b      An offset.
 * \note It is possible that vec_x = vec_y on input to overwrite x
 */
template <std::size_t dim>
constexpr void
t8_axb (const t8_vec<dim> &vec_x, t8_vec<dim> &vec_y, const double alpha, const double b)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), [alpha, b] (double v) { return alpha * v + b; });
}

/** Y = Y + alpha * X
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [in,out] vec_y On input, An N-dimensional vector.
 *                      On output set \a to vec_y + \a alpha * \a vec_x
 * \param [in]  alpha  A factor.
 */
template <std::size_t dim>
constexpr void
t8_axpy (const t8_vec<dim> &vec_x, t8_vec<dim> &vec_y, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), vec_y.begin (),
                  [alpha] (double x, double y) { return y + alpha * x; });
}

/** Z = Y + alpha * X
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [in]  vec_y  An N-dimensional vector.
 * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
 * \param [in]  alpha  A factor for the multiplication of \a vec_x.
 */
template <std::size_t dim>
constexpr void
t8_axpyz (const t8_vec<dim> &vec_x, const t8_vec<dim> &vec_y, t8_vec<dim> &vec_z, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), vec_z.begin (),
                  [alpha] (double x, double y) { return y + alpha * x; });
}

/** Dot product of X and Y.
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [in]  vec_y  An N-dimensional vector.
 * \return             The dot product \a vec_x * \a vec_y
 */
template <std::size_t dim>
constexpr double
t8_dot (const t8_vec<dim> &vec_x, const t8_vec<dim> &vec_y)
{
  return std::inner_product (vec_x.begin (), vec_x.end (), vec_y.begin (), 0.0);
}
/** Cross product of X and Y
 * \param [in]  vec_x  A 2D vector.
 * \param [in]  vec_y  A 2D vector.
 * \return             The cross product of \a vec_x and \a vec_y.
 */
static inline double
t8_cross_2D (const t8_vec<2> &vec_x, const t8_vec<2> &vec_y)
{
  return vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Cross product of X and Y
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
static inline void
t8_cross_3D (const t8_3D_vec &vec_x, const t8_3D_vec &vec_y, t8_3D_vec &cross)
{
  cross[0] = vec_x[1] * vec_y[2] - vec_x[2] * vec_y[1];
  cross[1] = vec_x[2] * vec_y[0] - vec_x[0] * vec_y[2];
  cross[2] = vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Compute the difference of two vectors.
 * \param [in]  vec_x  An N-dimensional vector.
 * \param [in]  vec_y  An N-dimensional vector.
 * \param [out] diff   On output, the difference of \a vec_x and \a vec_y.
 */
template <std::size_t dim>
constexpr void
t8_diff (const t8_vec<dim> &vec_x, const t8_vec<dim> &vec_y, t8_vec<dim> &diff)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), diff.begin (), std::minus<double> ());
}

/**
 * Check the equality of two vectors or points elementwise 
 * 
 * \param[in] vec_x 
 * \param[in] vec_y 
 * \param[in] tol 
 * \return true, if the vectors are equal up to \a tol 
 */
template <typename T>
constexpr bool
t8_eq (const T &vec_x, const T &vec_y, const double tol)
{
  return std::equal (vec_x.begin (), vec_x.end (), vec_y.begin (),
                     [tol] (double x, double y) { return std::fabs (x - y) <= tol; });
}

/** Rescale a vector to a new length.
 * \param [in,out] vec  An N-dimensional vector.
 * \param [in]  new_length  New length of the vector.
 */
template <std::size_t dim>
static inline void
t8_rescale (t8_vec<dim> &vec, const double new_length)
{
  t8_normalize (vec);
  t8_ax (vec, new_length);
}

/** Compute the normal of a triangle given by its three vertices.
 * \param [in]  p1  A 3D vector.
 * \param [in]  p2  A 3D vector.
 * \param [in]  p3  A 3D vector.
 * \param [out] normal vector of the triangle. (Not necessarily of length 1!)
 */
static inline void
t8_normal_of_tri (const t8_3D_vec &p1, const t8_3D_vec &p2, const t8_3D_vec &p3, t8_3D_vec &normal)
{
  t8_3D_vec a;
  t8_3D_vec b;
  std::transform (p2.begin (), p2.end (), p1.begin (), a.begin (), std::minus<double> ());
  std::transform (p3.begin (), p3.end (), p1.begin (), b.begin (), std::minus<double> ());
  t8_cross_3D (a, b, normal);
}

/** Compute an orthogonal coordinate system from a given vector.
 * \param [in]   v1 3D vector.
 * \param [out]  v2 3D vector.
 * \param [out]  v3 3D vector.
 */
static inline void
t8_orthogonal_tripod (const t8_3D_vec &v1, t8_3D_vec &v2, t8_3D_vec &v3)
{
  v2[0] = v1[1];
  v2[1] = v1[2];
  v2[2] = -v1[0];

  t8_axpy<3> (v1, v2, -t8_dot<3> (v1, v2));
  t8_cross_3D (v1, v2, v3);

  t8_normalize<3> (v2);
  t8_normalize<3> (v3);
}

#endif /* !T8_VEC_HXX */
