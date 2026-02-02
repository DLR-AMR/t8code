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

#include <algorithm>
#include <numeric>
#include <span>
#include <array>
#include <type_traits>
#include <iterator>
#include <concepts>
#include <cmath>
#include <functional>

template <std::size_t TDim, typename TType = double>
using t8_vec = std::array<TType, TDim>;

/** Type alias for a 2D point.
 */
typedef std::array<double, 2> t8_2D_point;

/** Type alias for a 3D point.
 */
typedef std::array<double, 3> t8_3D_point;

/** Type alias for a 2D vector.
 */
typedef t8_vec<2> t8_2D_vec;

/** Type alias for a 3D vector.
 */
typedef t8_vec<3> t8_3D_vec;

/** Type alias for a non-owning 3D vec view.
 * \tparam TType The type (const and so on)
 */
template <std::size_t TDim, typename TType = double>
using t8_vec_view = std::span<TType, TDim>;

/** Convenience function to create a vector view from a raw array.
 *
 * \tparam TDim            The dimension of the array.
 * \tparam TType           The type (const and so on)
 * \param [in] ptr         The pointer to the array.
 * \return                 The view.
 */
template <std::size_t TDim, typename TType = double>
constexpr auto
make_t8_vec_view (TType *ptr) noexcept
{
  using view_type = t8_vec_view<TDim, TType>;
  return view_type (ptr, TDim);
}

/** Convenience function to create a 2D vector view from a raw array.
 *
 * \param [in] ptr         The pointer to the array.
 * \tparam TType           The type (const and so on)
 * \return                 The view.
 */
template <typename TType>
inline auto
make_t8_2D_vec_view (TType *ptr) noexcept
{
  return make_t8_vec_view<2> (ptr);
}

/** Convenience function to create a 3D vector view from a raw array.
 *
 * \param [in] ptr         The pointer to the array.
 * \tparam TType           The type (const and so on)
 * \return                 The view.
 */
template <typename TType>
inline auto
make_t8_3D_vec_view (TType *ptr) noexcept
{
  return make_t8_vec_view<3> (ptr);
}

/** Concept for container types with value type that is convertible into double.
*/
template <typename T>
concept T8ContainerdoubleType = requires (T t) {
  {
    std::begin (t)
  } -> std::input_iterator;
  {
    std::end (t)
  } -> std::input_iterator;
  typename T::value_type;
} && std::is_convertible_v<typename T::value_type, double>;

/** Concept for container types with any value type.
*/
template <typename T>
concept T8ContainerType = requires (T t) {
  {
    std::begin (t)
  } -> std::input_iterator;
  {
    std::end (t)
  } -> std::input_iterator;
};

/** Concept for container types of size 2 with value type that is convertible into double.
*/
template <typename T>
concept T8ContainerdoubleType2D = T8ContainerdoubleType<T> && requires (T t) {
  {
    std::size (t)
  } -> std::convertible_to<std::size_t>;
  requires std::size (t) == 2;
};

/** Concept for container types of size 3 with value type that is convertible into double.
*/
template <typename T>
concept T8ContainerdoubleType3D = T8ContainerdoubleType<T> && requires (T t) {
  {
    std::size (t)
  } -> std::convertible_to<std::size_t>;
  requires std::size (t) == 3;
};

/** Vector norm.
  * \param [in] vec  An N-dimensional vector.
  * \return          The norm of \a vec.
  */
template <T8ContainerdoubleType TVec>
static inline double
t8_norm (const TVec &vec)
{
  return std::sqrt (std::inner_product (vec.begin (), vec.end (), vec.begin (), 0.0));
}

/** Normalize a vector.
  * \param [in,out] vec  An N-dimensional vector.
  */
template <T8ContainerdoubleType TVec>
constexpr void
t8_normalize (TVec &vec)
{
  const double norm = t8_norm (vec);
  std::transform (vec.begin (), vec.end (), vec.begin (), [norm] (double v) { return v / norm; });
}

/** Copy a dimensional object.
  * \param [in]  src  The source.
  * \param [out] dest The destination.
  */
template <T8ContainerType TVec1, T8ContainerType TVec2>
constexpr void
t8_copy (const TVec1 &src, TVec2 &dest)
{
  std::copy (src.begin (), src.end (), dest.begin ());
}

/** Euclidean distance of X and Y.
  * \param [in]  point_x  An N-dimensional point.
  * \param [in]  point_y  An N-dimensional point.
  * \return             The euclidean distance.
  *                     Equivalent to norm (X-Y).
  */
template <T8ContainerdoubleType TPointX, T8ContainerdoubleType TPointY>
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
template <T8ContainerdoubleType TVec>
constexpr void
t8_ax (TVec &vec_x, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_x.begin (), [alpha] (double v) { return v * alpha; });
}

/** Compute Y = alpha * X
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [out] vec_y  On output set to \a alpha * \a vec_x.
  * \param [in]  alpha  A factor.
  */
template <T8ContainerdoubleType TVecX, T8ContainerdoubleType TVecY>
constexpr void
t8_axy (const TVecX &vec_x, TVecY &vec_y, const double alpha)
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
template <T8ContainerdoubleType TVecX, T8ContainerdoubleType TVecY>
constexpr void
t8_axb (const TVecX &vec_x, TVecY &vec_y, const double alpha, const double b)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), [alpha, b] (double v) { return alpha * v + b; });
}

/** Y = Y + alpha * X
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in,out] vec_y On input, An N-dimensional vector.
  *                      On output set \a to vec_y + \a alpha * \a vec_x
  * \param [in]  alpha  A factor.
  */
template <T8ContainerdoubleType TVecX, T8ContainerdoubleType TVecY>
constexpr void
t8_axpy (const TVecX &vec_x, TVecY &vec_y, const double alpha)
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
template <T8ContainerdoubleType TVecX, T8ContainerdoubleType TVecY, T8ContainerdoubleType TVecZ>
constexpr void
t8_axpyz (const TVecX &vec_x, const TVecY &vec_y, TVecZ &vec_z, const double alpha)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), vec_z.begin (),
                  [alpha] (double x, double y) { return y + alpha * x; });
}

/** Dot product of X and Y.
  * \param [in]  vec_x  An N-dimensional vector.
  * \param [in]  vec_y  An N-dimensional vector.
  * \return             The dot product \a vec_x * \a vec_y
  */
template <T8ContainerdoubleType TVecX, T8ContainerdoubleType TVecY>
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
template <T8ContainerdoubleType2D TVecX, T8ContainerdoubleType2D TVecY>
static inline double
t8_cross_2D (const TVecX &vec_x, const TVecY &vec_y)
{
  return vec_x[0] * vec_y[1] - vec_x[1] * vec_y[0];
}

/** Cross product of X and Y
  * \param [in]  vec_x  A 3D vector.
  * \param [in]  vec_y  A 3D vector.
  * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
  */
template <T8ContainerdoubleType3D TVecX, T8ContainerdoubleType3D TVecY, T8ContainerdoubleType3D TVecCross>
static inline void
t8_cross_3D (const TVecX &vec_x, const TVecY &vec_y, TVecCross &cross)
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
template <T8ContainerdoubleType3D TVecX, T8ContainerdoubleType3D TVecY, T8ContainerdoubleType3D TVecDiff>
constexpr void
t8_diff (const TVecX &vec_x, const TVecY &vec_y, TVecDiff &diff)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), diff.begin (), std::minus<double> ());
}

/**
  * Check the equality of two dimensional objects.
  * \param[in] x Container 1 to compare.
  * \param[in] y Container 2 that should be compared to \a x for equality given tolerance \a tol.
  * \param[in] tol Tolerance.
  * \return true, if the objects are equal up to \a tol.
  */
template <T8ContainerdoubleType TDimensionalX, T8ContainerdoubleType TDimensionalY>
constexpr bool
t8_eq (const TDimensionalX &x, const TDimensionalY &y, const double tol)
{
  return std::equal (x.begin (), x.end (), y.begin (),
                     [tol] (double x_val, double y_val) { return std::fabs (x_val - y_val) <= tol; });
}

/** Rescale a vector to a new length.
  * \param [in,out] vec  An N-dimensional vector.
  * \param [in]  new_length  New length of the vector.
  */
template <T8ContainerdoubleType TVec>
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
  * \param [out] normal vector of the triangle. (Not necessarily of length 1!)d
  */
template <T8ContainerdoubleType3D TVecP1, T8ContainerdoubleType3D TVecP2, T8ContainerdoubleType3D TVecP3,
          T8ContainerdoubleType3D TVecNormal>
static inline void
t8_normal_of_tri (const TVecP1 &p1, const TVecP2 &p2, const TVecP3 &p3, TVecNormal &normal)
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
template <T8ContainerdoubleType3D TVecV1, T8ContainerdoubleType3D TVecV2, T8ContainerdoubleType3D TVecV3>
static inline void
t8_orthogonal_tripod (const TVecV1 &v1, TVecV2 &v2, TVecV3 &v3)
{
  v2[0] = v1[1];
  v2[1] = v1[2];
  v2[2] = -v1[0];

  t8_axpy (v1, v2, -t8_dot (v1, v2));
  t8_cross_3D (v1, v2, v3);

  t8_normalize (v2);
  t8_normalize (v3);
}

#endif /* !T8_VEC_HXX */
