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

#include <algorithm>
#include <numeric>
#include <span>

template <std::size_t TDim>
struct t8_vec_tag
{
  static constexpr std::size_t dim = TDim;
};

/**
  * Type alias for a vector in N-dimensional space.
  * \tparam TDim Dimension of the vector.
  */
template <std::size_t TDim>
using t8_vec = T8Type<std::array<double, TDim>, t8_vec_tag<TDim>, EqualityComparable, Swapable, RandomAccessible>;

/**
  * Type alias for a non-owning vector view in N-dimensional space.
  * \tparam TDim     Dimension of the vector.
  * \tparam TType    The type (const and so on)
  */
template <std::size_t TDim, typename TType = double>
using t8_vec_view = T8Type<std::span<TType, TDim>, t8_vec_tag<TDim>, EqualityComparable, Swapable, RandomAccessible>;

/** Type alias for a 2D vector.
 */
using t8_2D_vec = t8_vec<2>;

/** Type alias for a non-owning 2D vec view.
 * \tparam TType The type (const and so on)
 */
template <typename TType = double>
using t8_2D_vec_view = t8_vec_view<2, TType>;

/** Type alias for a 3D vector.
 */
using t8_3D_vec = t8_vec<3>;

/** Type alias for a non-owning 3D vec view.
 * \tparam TType The type (const and so on)
 */
template <typename TType = double>
using t8_3D_vec_view = t8_vec_view<3, TType>;

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
  return view_type (std::span<TType, TDim> (ptr, TDim));
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

/* -----------------------Concepts for t8_vec----------------------- */

/** Concept that checks whether a type is a strong type of type t8_vec<N>. N can be either fixed or left open.
  * \tparam TType          The type to check.
  * \tparam TExpectedDim   Optional dimensional restriction (default = wildcard).
  */
template <typename T, std::size_t Expected = static_cast<std::size_t> (-1)>
concept T8VecType = requires { typename std::remove_cvref_t<T>::tag; } && requires {
  {
    std::remove_cvref_t<T>::tag::dim
  } -> std::convertible_to<std::size_t>;
} && (Expected == static_cast<std::size_t> (-1) || std::remove_cvref_t<T>::tag::dim == Expected);

/** Concept that checks whether a type is a container of elements of type t8_vec<N>. N can be either fixed or left open.
  * \tparam TType     The type to check.
  * \tparam TDim   Optional dimensional restriction (default = wildcard).
  */
template <typename TType, std::size_t TExpectedDim = static_cast<std::size_t> (-1)>
concept T8VecContainerType = std::ranges::range<std::remove_cvref_t<TType>>
                             && T8VecType<std::ranges::range_value_t<std::remove_cvref_t<TType>>, TExpectedDim>;

/* -----------------------End concepts for t8_vec----------------------- */

template <std::size_t TDim>
struct t8_point_tag
{
  static constexpr std::size_t dim = TDim;
};

/**
  * Type alias for a point in N-dimensional space.
  * \tparam TDim Dimension of the point.
  */
template <std::size_t TDim>
using t8_point = T8Type<std::array<double, TDim>, t8_point_tag<TDim>, EqualityComparable, Swapable, RandomAccessible>;

/**
  * Type alias for a non-owning point view in N-dimensional space.
  * \tparam TDim     Dimension of the point.
  * \tparam TType    The type (const and so on)
  */
template <std::size_t TDim, typename TType = double>
using t8_point_view
  = T8Type<std::span<TType, TDim>, t8_point_tag<TDim>, EqualityComparable, Swapable, RandomAccessible>;

/** Type alias for a 2D point.
 */
using t8_2D_point = t8_point<2>;

/** Type alias for a non-owning 2D point view.
 * \tparam TType The type (const and so on)
 */
template <typename TType = double>
using t8_2D_point_view = t8_point_view<2, TType>;

/** Type alias for a 3D point.
 */
using t8_3D_point = t8_point<3>;

/** Type alias for a non-owning 3D point view.
 * \tparam TType The type (const and so on)
 */
template <typename TType = double>
using t8_3D_point_view = t8_point_view<3, TType>;

/** Convenience function to create a point view from a raw array.
 *
 * \tparam TDim            The dimension of the array.
 * \tparam TType           The type (const and so on)
 * \param [in] ptr         The pointer to the array.
 * \return                 The view.
 */
template <std::size_t TDim, typename TType = double>
constexpr auto
make_t8_point_view (TType *ptr) noexcept
{
  using view_type = t8_point_view<TDim, TType>;
  return view_type (std::span<TType, TDim> (ptr, TDim));
}

/** Convenience function to create a 2D point view from a raw array.
 *
 * \param [in] ptr         The pointer to the array.
 * \tparam TType           The type (const and so on)
 * \return                 The view.
 */
template <typename TType>
inline auto
make_t8_2D_point_view (TType *ptr) noexcept
{
  return make_t8_point_view<2> (ptr);
}

/** Convenience function to create a 3D point view from a raw array.
 *
 * \param [in] ptr         The pointer to the array.
 * \tparam TType           The type (const and so on)
 * \return                 The view.
 */
template <typename TType>
inline auto
make_t8_3D_point_view (TType *ptr) noexcept
{
  return make_t8_point_view<3> (ptr);
}

/* -----------------------Concepts for t8_point----------------------- */
/** Concept that checks whether a type is a strong type of type t8_point<N>.
  * N can be either fixed or left open.
  * \tparam TType          The type to check.
  * \tparam TExpectedDim   Optional dimensional restriction (default = wildcard).
  */
template <typename T, std::size_t Expected = static_cast<std::size_t> (-1)>
concept T8PointType = requires { typename std::remove_cvref_t<T>::tag; } && requires {
  {
    std::remove_cvref_t<T>::tag::dim
  } -> std::convertible_to<std::size_t>;
} && (Expected == static_cast<std::size_t> (-1) || std::remove_cvref_t<T>::tag::dim == Expected);

/** Concept that checks whether a type is a container of elements of type t8_point<N>.
   * N can be either fixed or left open.
   * \tparam TType          The type to check.
   * \tparam TExpectedDim   Optional dimensional restriction (default = wildcard).
   */
template <typename TType, std::size_t TExpectedDim = static_cast<std::size_t> (-1)>
concept T8PointContainerType = std::ranges::range<std::remove_cvref_t<TType>>
                               && T8PointType<std::ranges::range_value_t<std::remove_cvref_t<TType>>, TExpectedDim>;

/* -----------------------End concepts for t8_point----------------------- */

/** Concept that checks if a type is either a point or vec
  * \tparam TType     The type to check.
  * \tparam TDim   Optional dimensional restriction (default = wildcard).
  */
template <typename TType, std::size_t TDim = static_cast<std::size_t> (-1)>
concept T8DimensionalType = T8PointType<TType, TDim> || T8VecType<TType, TDim>;

/**
  * Dimension of a dimensional type.
  * \tparam TDimensional     The dimensional type.
  * \return                  The dimension.
  */
template <T8DimensionalType TDimensional>
constexpr std::size_t dim_v = std::remove_cvref_t<TDimensional>::tag_type::dim;

/** Vector norm.
  * \param [in] vec  An N-dimensional vector.
  * \return          The norm of \a vec.
  */
template <T8VecType TVec>
static inline double
t8_norm (const TVec &vec)
{
  return std::sqrt (std::inner_product (vec.begin (), vec.end (), vec.begin (), 0.0));
}

/** Normalize a vector.
  * \param [in,out] vec  An N-dimensional vector.
  */
template <T8VecType TVec>
constexpr void
t8_normalize (TVec &vec)
{
  const double norm = t8_norm (vec);
  std::transform (vec.begin (), vec.end (), vec.begin (), [norm] (double v) { return v / norm; });
}

/**
  * Copy a dimensional object.
  * \param [in]  src  The source.
  * \param [out] dest The destination.
  */
template <T8DimensionalType TDimensional1, T8DimensionalType TDimensional2>
constexpr void
t8_copy (const TDimensional1 &src, TDimensional2 &dest)
{
  std::copy (src.begin (), src.end (), dest.begin ());
}

/** Euclidean distance of X and Y.
  * \param [in]  point_x  An N-dimensional point.
  * \param [in]  point_y  An N-dimensional point.
  * \return             The euclidean distance.
  *                     Equivalent to norm (X-Y).
  */
template <T8PointType TPointX, T8PointType TPointY>
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
template <T8VecType TVec>
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
template <T8VecType TVecX, T8VecType TVecY>
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
template <T8VecType TVecX, T8VecType TVecY>
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
template <T8VecType TVecX, T8VecType TVecY>
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
template <T8VecType TVecX, T8VecType TVecY, T8VecType TVecZ>
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
template <T8VecType TVecX, T8VecType TVecY>
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
template <T8VecType<2> TVecX, T8VecType<2> TVecY>
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
template <T8VecType<3> TVecX, T8VecType<3> TVecY, T8VecType<3> TVecCross>
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

template <T8VecType<3> TVecX, T8VecType<3> TVecY, T8VecType<3> TVecDiff>
constexpr void
t8_diff (const TVecX &vec_x, const TVecY &vec_y, TVecDiff &diff)
{
  std::transform (vec_x.begin (), vec_x.end (), vec_y.begin (), diff.begin (), std::minus<double> ());
}

/**
  * Check the equality of two dimensional objects
  * \param[in] x
  * \param[in] y
  * \param[in] tol
  * \return true, if the objects are equal up to \a tol
  */
template <T8DimensionalType TDimensionalX, T8DimensionalType TDimensionalY>
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
template <T8VecType TVec>
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

template <T8VecType<3> TVecP1, T8VecType<3> TVecP2, T8VecType<3> TVecP3, T8VecType<3> TVecNormal>
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
template <T8VecType<3> TVecV1, T8VecType<3> TVecV2, T8VecType<3> TVecV3>
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
