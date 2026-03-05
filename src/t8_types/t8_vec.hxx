/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
 * We define routines to handle vectors.
 */

#ifndef T8_VEC_HXX
#define T8_VEC_HXX

#include <t8.h>
#include <array>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <concepts>
#include <cmath>
#include <functional>

/** Type alias for a vector in N-dimensional space.
 * \tparam TDim Dimension of the vector.
 * \tparam TType Value type of the vector.
 */
template <std::size_t TDim, typename TType = double>
using t8_vec = std::array<TType, TDim>;

/** Type alias for a 2D vector.
 */
using t8_2D_vec = t8_vec<2>;

/** Type alias for a 3D vector.
 */
using t8_3D_vec = t8_vec<3>;

/** Concept for input ranges with elements convertible to double.
 * \tparam TType Container type to check.
 */
template <typename TType>
concept T8InputRange
  = std::ranges::input_range<TType> && std::convertible_to<std::ranges::range_value_t<TType>, double>;

/** Concept for random access ranges with elements convertible to double.
 * \note T8RandomAccessRange also satisfies the condition of T8InputRange but not vice versa.
 * \tparam TType Container type to check.
 */
template <typename TType>
concept T8RandomAccessRange
  = std::ranges::random_access_range<TType> && std::convertible_to<std::ranges::range_value_t<TType>, double>;

/** Vector norm.
  * \param [in] vec  An N-dimensional vector.
  * \return          The norm of \a vec.
  */
template <T8InputRange TVec>
static inline double
t8_norm (const TVec &vec)
{
  return std::sqrt (std::inner_product (vec.begin (), vec.end (), vec.begin (), 0.0));
}

/** Normalize a vector.
  * \param [in,out] vec  An N-dimensional vector.
  */
template <T8InputRange TVec>
constexpr void
t8_normalize (TVec &vec)
{
  const double norm = t8_norm (vec);
  T8_ASSERT (norm != 0);
  std::ranges::transform (vec, vec.begin (), [norm] (double v) { return v / norm; });
}

/** Copy a dimensional object.
  * \param [in]  src  The source.
  * \param [out] dest The destination.
  */
template <T8InputRange TVec1, T8InputRange TVec2>
constexpr void
t8_copy (const TVec1 &src, TVec2 &dest)
{
  std::ranges::copy (src, dest.begin ());
}

/** Euclidean distance of X and Y.
  * \param [in]  point_x  An N-dimensional point.
  * \param [in]  point_y  An N-dimensional point.
  * \return             The euclidean distance.
  *                     Equivalent to norm (X-Y).
  */
template <T8InputRange TPointX, T8InputRange TPointY>
constexpr double
t8_dist (const TPointX &point_x, const TPointY &point_y)
{
  double dist = std::inner_product (point_x.begin (), point_x.end (), point_y.begin (), 0.0, std::plus<double> (),
                                    [] (double x, double y) { return (x - y) * (x - y); });
  return std::sqrt (dist);
}

/** Compute X = alpha * X
  * \param [in,out] vec_x  An N-dimensional vector. On output set to \a alpha * \a vec_x.
  * \param [in]     alpha  A factor.
  */
template <T8InputRange TVec>
constexpr void
t8_ax (TVec &vec_x, const double alpha)
{
  std::ranges::transform (vec_x, vec_x.begin (), [alpha] (double v) { return v * alpha; });
}

/** Compute Y = alpha * X
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [out] vec_y  On output set to \a alpha * \a vec_x.
  * \param [in]  alpha  A factor.
  */
template <T8InputRange TVecX, T8InputRange TVecY>
constexpr void
t8_axy (const TVecX &vec_x, TVecY &vec_y, const double alpha)
{
  std::ranges::transform (vec_x, vec_y.begin (), [alpha] (double v) { return v * alpha; });
}

/** Y = alpha * X + b
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [out] vec_y  On input, An N-dimensional vector.
  *                     On output set to \a alpha * \a vec_x + \a b.
  * \param [in]  alpha  A factor.
  * \param [in]  b      An offset.
  */
template <T8InputRange TVecX, T8InputRange TVecY>
constexpr void
t8_axb (const TVecX &vec_x, TVecY &vec_y, const double alpha, const double b)
{
  std::ranges::transform (vec_x, vec_y.begin (), [alpha, b] (double v) { return alpha * v + b; });
}

/** Y = Y + alpha * X
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in,out] vec_y On input, An N-dimensional vector.
  *                      On output set \a to vec_y + \a alpha * \a vec_x
  * \param [in]  alpha  A factor.
  */
template <T8InputRange TVecX, T8InputRange TVecY>
constexpr void
t8_axpy (const TVecX &vec_x, TVecY &vec_y, const double alpha)
{
  std::ranges::transform (vec_x, vec_y, vec_y.begin (), [alpha] (double x, double y) { return y + alpha * x; });
}

/** Z = Y + alpha * X
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in]  vec_y  An N-dimensional vector.
  * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
  * \param [in]  alpha  A factor for the multiplication of \a vec_x.
  */
template <T8InputRange TVecX, T8InputRange TVecY, T8InputRange TVecZ>
constexpr void
t8_axpyz (const TVecX &vec_x, const TVecY &vec_y, TVecZ &vec_z, const double alpha)
{
  std::ranges::transform (vec_x, vec_y, vec_z.begin (), [alpha] (double x, double y) { return y + alpha * x; });
}

/** Dot product of X and Y.
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in]  vec_y  An N-dimensional vector.
  * \return             The dot product \a vec_x * \a vec_y
  */
template <T8InputRange TVecX, T8InputRange TVecY>
constexpr double
t8_dot (const TVecX &vec_x, const TVecY &vec_y)
{
  return std::inner_product (vec_x.begin (), vec_x.end (), vec_y.begin (), 0.0);
}

/** Cross product of X and Y
  * \param [in]  vec_x  A 2D vector.
  * \param [in]  vec_y  A 2D vector.
  * \return             The cross product of \a vec_x and \a vec_y.
  */
template <T8RandomAccessRange TVecX, T8RandomAccessRange TVecY>
static inline double
t8_cross_2D (const TVecX &vec_x, const TVecY &vec_y)
{
  T8_ASSERT ((std::ranges::distance (vec_x) >= 2) && (std::ranges::distance (vec_y) >= 2));
  return vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Cross product of X and Y
  * \param [in]  vec_x  A 3D vector.
  * \param [in]  vec_y  A 3D vector.
  * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
  */
template <T8RandomAccessRange TVecX, T8RandomAccessRange TVecY, T8RandomAccessRange TVecCross>
static inline void
t8_cross_3D (const TVecX &vec_x, const TVecY &vec_y, TVecCross &cross)
{
  T8_ASSERT ((std::ranges::distance (vec_x) >= 3) && (std::ranges::distance (vec_y) >= 3)
             && (std::ranges::distance (cross) >= 3));
  cross[0] = vec_x[1] * vec_y[2] - vec_x[2] * vec_y[1];
  cross[1] = vec_x[2] * vec_y[0] - vec_x[0] * vec_y[2];
  cross[2] = vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Compute the difference of two vectors.
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in]  vec_y  An N-dimensional vector.
  * \param [out] diff   On output, the difference of \a vec_x and \a vec_y.
  */
template <T8InputRange TVecX, T8InputRange TVecY, T8InputRange TVecDiff>
constexpr void
t8_diff (const TVecX &vec_x, const TVecY &vec_y, TVecDiff &diff)
{
  T8_ASSERT (std::ranges::distance (vec_x) == std::ranges::distance (vec_y));
  std::ranges::transform (vec_x, vec_y, diff.begin (), std::minus {});
}

/**
  * Check the equality of two dimensional objects.
  * \param[in] x Container 1 to compare.
  * \param[in] y Container 2 that should be compared to \a x for equality given tolerance \a tol.
  * \param[in] tol Tolerance.
  * \return true, if the objects are equal up to \a tol.
  */
template <T8InputRange TDimensionalX, T8InputRange TDimensionalY>
constexpr bool
t8_eq (const TDimensionalX &x, const TDimensionalY &y, const double tol)
{
  return std::ranges::equal (x, y, [tol] (double x_val, double y_val) { return std::fabs (x_val - y_val) <= tol; });
}

/** Rescale a vector to a new length.
  * \param [in,out] vec  An N-dimensional vector.
  * \param [in]  new_length  New length of the vector.
  */
template <T8InputRange TVec>
static inline void
t8_rescale (TVec &vec, const double new_length)
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
template <T8InputRange TVecP1, T8InputRange TVecP2, T8InputRange TVecP3, T8InputRange TVecNormal>
static inline void
t8_normal_of_tri (const TVecP1 &p1, const TVecP2 &p2, const TVecP3 &p3, TVecNormal &normal)
{
  T8_ASSERT ((std::ranges::distance (p1) >= 3) && (std::ranges::distance (p2) >= 3)
             && (std::ranges::distance (p3) >= 3));

  t8_3D_vec a;
  t8_3D_vec b;

  std::ranges::transform (p2, p1, a.begin (), std::minus {});
  std::ranges::transform (p3, p1, b.begin (), std::minus {});

  t8_cross_3D (a, b, normal);
}

/** Compute an orthogonal coordinate system from a given vector.
  * \param [in]   v1 3D vector.
  * \param [out]  v2 3D vector.
  * \param [out]  v3 3D vector.
  */
template <T8RandomAccessRange TVecV1, T8RandomAccessRange TVecV2, T8RandomAccessRange TVecV3>
static inline void
t8_orthogonal_tripod (const TVecV1 &v1, TVecV2 &v2, TVecV3 &v3)
{
  T8_ASSERT ((std::ranges::distance (v1) >= 3) && (std::ranges::distance (v2) >= 3)
             && (std::ranges::distance (v3) >= 3));
  v2[0] = v1[1];
  v2[1] = v1[2];
  v2[2] = -v1[0];

  t8_axpy (v1, v2, -t8_dot (v1, v2));
  t8_cross_3D (v1, v2, v3);

  t8_normalize (v2);
  t8_normalize (v3);
}

#endif /* !T8_VEC_HXX */
