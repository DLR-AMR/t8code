/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/* In this file we collect tests for the routines in t8_vec.hxx */

#include <gtest/gtest.h>
#include <t8_types/t8_vec.hxx>
#include <t8_types/t8_vec.h>
#include <t8_gtest_custom_assertion.hxx>
#include <t8_gtest_memory_macros.hxx>
#include <t8_helper_functions/t8_unrolled_for.hxx>

#include <random>

/* test the t8_norm function */
TEST (t8_gtest_vec, norm)
{
  const t8_3D_vec zero ({ 0, 0, 0 });
  const t8_3D_vec onetwothree ({ 1, 2, 3 });
  const t8_3D_vec arbitrary ({ -.05, 3.14159, 42 });

  const double normonetwothree = sqrt (1 + 4 + 9);
  const double normarbitrary = 42.117360883;

  EXPECT_EQ (t8_norm (zero), 0);
  EXPECT_NEAR (t8_norm (onetwothree), normonetwothree, T8_PRECISION_SQRT_EPS);
  EXPECT_NEAR (t8_norm (arbitrary), normarbitrary, T8_PRECISION_SQRT_EPS);
}

/* test the t8_dist function */
TEST (t8_gtest_vec, dist)
{
  const t8_3D_point zero ({ 0, 0, 0 });
  const t8_3D_point onetwothree ({ 1, 2, 3 });
  const t8_3D_point arbitrary ({ -.05, 3.14159, 42 });
  const double distzeroonetwothree = sqrt (1 + 4 + 9);
  const double distarbitraryonetwothree = 39.030830477;
  EXPECT_VEC_EQ (zero, zero, T8_PRECISION_SQRT_EPS);
  EXPECT_VEC_EQ (onetwothree, onetwothree, T8_PRECISION_SQRT_EPS);
  EXPECT_NEAR (t8_dist (onetwothree, zero), distzeroonetwothree, T8_PRECISION_SQRT_EPS);
  EXPECT_NEAR (t8_dist (zero, onetwothree), distzeroonetwothree, T8_PRECISION_SQRT_EPS);
  EXPECT_NEAR (t8_dist (arbitrary, onetwothree), distarbitraryonetwothree, T8_PRECISION_SQRT_EPS);
  EXPECT_NEAR (t8_dist (onetwothree, arbitrary), distarbitraryonetwothree, T8_PRECISION_SQRT_EPS);
}

/* test the t8_ax function */
TEST (t8_gtest_vec, ax)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  t8_3D_vec zero ({ 0, 0, 0 });
  t8_3D_vec onetwothree ({ 1, 2, 3 });
  t8_3D_vec arbitrary ({ -.05, 3.14159, 42 });
  const double alpha = 5.1234;

  /* Compute Y = alpha * Y */
  t8_ax (zero, alpha);
  t8_ax (onetwothree, alpha);
  t8_ax (arbitrary, alpha);

  /* Check results */
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (zero[i], alpha * czero[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (onetwothree[i], alpha * conetwothree[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (arbitrary[i], alpha * carbitrary[i], T8_PRECISION_SQRT_EPS);
  }
}

/* test the t8_axy function */
TEST (t8_gtest_vec, axy)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  t8_3D_vec zero;
  t8_3D_vec onetwothree;
  t8_3D_vec arbitrary;
  const double alpha = 5.1234;

  /* Compute Y = alpha * X */
  t8_axy (czero, zero, alpha);
  t8_axy (conetwothree, onetwothree, alpha);
  t8_axy (carbitrary, arbitrary, alpha);

  /* Check results */
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (zero[i], alpha * czero[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (onetwothree[i], alpha * conetwothree[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (arbitrary[i], alpha * carbitrary[i], T8_PRECISION_SQRT_EPS);
  }
}

/* test the t8_axb function */
TEST (t8_gtest_vec, axb)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  const double b = 2.71828;
  t8_3D_vec zero;
  t8_3D_vec onetwothree;
  t8_3D_vec arbitrary;
  const double alpha = 5.1234;

  /* Compute Y = alpha * Y + b */
  t8_axb (czero, zero, alpha, b);
  t8_axb (conetwothree, onetwothree, alpha, b);
  t8_axb (carbitrary, arbitrary, alpha, b);

  /* Check results */
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (zero[i], alpha * czero[i] + b, T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (onetwothree[i], alpha * conetwothree[i] + b, T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (arbitrary[i], alpha * carbitrary[i] + b, T8_PRECISION_SQRT_EPS);
  }
}

/* test the t8_axpy function */
TEST (t8_gtest_vec, axpy)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  const t8_3D_vec init ({ 3, 2.71828, -4.1 });
  /* The next three vecs must be initialized with the values of init */
  t8_3D_vec zero ({ 3, 2.71828, -4.1 });
  t8_3D_vec onetwothree ({ 3, 2.71828, -4.1 });
  t8_3D_vec arbitrary ({ 3, 2.71828, -4.1 });
  const double alpha = 5.1234;

  /* Compute Y = Y + alpha * Y */
  t8_axpy (czero, zero, alpha);
  t8_axpy (conetwothree, onetwothree, alpha);
  t8_axpy (carbitrary, arbitrary, alpha);

  /* Check results */
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (zero[i], init[i] + alpha * czero[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (onetwothree[i], init[i] + alpha * conetwothree[i], T8_PRECISION_SQRT_EPS);
    EXPECT_NEAR (arbitrary[i], init[i] + alpha * carbitrary[i], T8_PRECISION_SQRT_EPS);
  }
}

/* test the t8_axpyz function */
TEST (t8_gtest_vec, axpyz)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  const t8_3D_vec init ({ 3, 2.71828, -4.1 });
  t8_3D_vec Z;
  const double alpha = 5.1234;

  /* Z = init + alpha * zero */
  t8_axpyz (czero, init, Z, alpha);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * czero[i], T8_PRECISION_SQRT_EPS);
  }
  /* Z = init + alpha * conetwothree */
  t8_axpyz (conetwothree, init, Z, alpha);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * conetwothree[i], T8_PRECISION_SQRT_EPS);
  }
  /* Z = init + alpha * carbitrary */
  t8_axpyz (carbitrary, init, Z, alpha);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * carbitrary[i], T8_PRECISION_SQRT_EPS);
  }
}

/* test the t8_dot function */
TEST (t8_gtest_vec, dot)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec conetwothree ({ 1, 2, 3 });
  const t8_3D_vec carbitrary ({ -.05, 3.14159, 42 });
  double result;

  /* Dot product with 0 is 0 */
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ (t8_dot (czero, czero), 0);
    EXPECT_EQ (t8_dot (czero, conetwothree), 0);
    EXPECT_EQ (t8_dot (czero, carbitrary), 0);
  }

  /* Test dot product of the remaining two vecs */
  result = 0;
  for (int i = 0; i < 3; ++i) {
    result += conetwothree[i] * carbitrary[i];
  }
  EXPECT_NEAR (t8_dot (conetwothree, carbitrary), result, T8_PRECISION_SQRT_EPS);

  /* For the dot-product of a vector with itself we use the square of its norm */
  result = t8_norm (conetwothree) * t8_norm (conetwothree);
  EXPECT_NEAR (t8_dot (conetwothree, conetwothree), result, T8_PRECISION_SQRT_EPS);

  result = t8_norm (carbitrary) * t8_norm (carbitrary);
  EXPECT_NEAR (t8_dot (carbitrary, carbitrary), result, T8_PRECISION_SQRT_EPS);
}

/* test the t8_cross_3D function */
TEST (t8_gtest_vec, cross_3D)
{
  const t8_3D_vec czero ({ 0, 0, 0 });
  const t8_3D_vec e1 ({ 1, 0, 0 });
  const t8_3D_vec e2 ({ 0, 1, 0 });
  const t8_3D_vec e3 ({ 0, 0, 1 });
  t8_3D_vec cross;

  /* cross product with 0 is 0 */
  t8_cross_3D (czero, czero, cross);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ (cross[i], 0);
  }
  t8_cross_3D (e1, czero, cross);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ (cross[i], 0);
  }
  t8_cross_3D (e2, czero, cross);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ (cross[i], 0);
  }

  /* e1 x e2 = e3 */
  t8_cross_3D (e1, e2, cross);
  EXPECT_VEC_EQ (cross, e3, T8_PRECISION_SQRT_EPS);

  /* e2 x e3 = e1 */
  t8_cross_3D (e2, e3, cross);
  EXPECT_VEC_EQ (cross, e1, T8_PRECISION_SQRT_EPS);
}

TEST (t8_gtest_vec, cross_2D)
{
  const t8_vec<2> zero ({ 0, 0 });
  const t8_vec<2> e1 ({ 1, 0 });   // Unit vector along x-axis
  const t8_vec<2> e2 ({ 0, 1 });   // Unit vector along y-axis
  const t8_vec<2> v1 ({ 3, 4 });   // Arbitrary vector
  const t8_vec<2> v2 ({ -4, 3 });  // Perpendicular to v1

  double cross;

  // Cross product with zero vector is 0
  cross = t8_cross_2D (zero, zero);
  EXPECT_EQ (cross, 0.0);

  cross = t8_cross_2D (e1, zero);
  EXPECT_EQ (cross, 0.0);

  cross = t8_cross_2D (e2, zero);
  EXPECT_EQ (cross, 0.0);

  // Cross product of e1 and e2 is +1 (counterclockwise)
  cross = t8_cross_2D (e1, e2);
  EXPECT_EQ (cross, 1.0);

  // Cross product of e2 and e1 is -1 (clockwise)
  cross = t8_cross_2D (e2, e1);
  EXPECT_EQ (cross, -1.0);

  // Cross product of v1 and v2
  cross = t8_cross_2D (v1, v2);
  EXPECT_EQ (cross, 3 * 3 - 4 * -4);  // 9 + 16 = 25
  EXPECT_EQ (cross, 25.0);

  // Cross product of v2 and v1 (reverse order)
  cross = t8_cross_2D (v2, v1);
  EXPECT_EQ (cross, -25.0);
}

TEST (t8_gtest_vec, check_less_or_equal)
{
  const t8_3D_vec one ({ 1.0, 1.0, 1.0 });
  const t8_3D_vec one_minus_eps (
    { 1.0 - T8_PRECISION_SQRT_EPS, 1.0 - T8_PRECISION_SQRT_EPS, 1.0 - T8_PRECISION_SQRT_EPS });

  EXPECT_VEC_EQ (one, one_minus_eps, T8_PRECISION_SQRT_EPS);
}

/**
 * Creates a vector of vec views for plain c vectors.
 * \tparam TDim                   The dimension of the vector.
 * \param [in, out] c_vectors     Pointer to the raw c vectors.
 * \param [in, out] num_vectors   Number of vectors.
 * \return                        A std::vector containing views to the c vectors.
 */
template <size_t TDim>
static inline std::vector<t8_vec_view<TDim, const double>>
t8_convert_array_to_vec_view (const double* c_vectors, const size_t num_vectors)
{
  std::vector<t8_vec_view<TDim, const double>> vec_views;
  vec_views.reserve (num_vectors);
  for (size_t ivec = 0; ivec < num_vectors; ++ivec)
    vec_views.emplace_back (make_t8_vec_view<TDim, const double> (c_vectors + ivec * TDim));
  T8_ASSERT (vec_views.size () == num_vectors);
  return vec_views;
}

/**
 * Create a vector of t8_vec as a copy of plain c vectors.
 * \tparam TDim                   The dimension of the vector.
 * \param [in, out] c_vectors     Pointer to the raw c vectors.
 * \param [in, out] num_vectors   Number of vectors.
 * \return                        A std::vector containing t8_vec objects.
 */
template <size_t TDim>
static inline std::vector<t8_vec<TDim>>
t8_convert_array_to_vec (const double* c_vectors, const size_t num_vectors)
{
  std::vector<t8_vec<TDim>> vecs (num_vectors);
  for (size_t ivec = 0; ivec < num_vectors; ++ivec)
    std::copy_n (c_vectors + ivec * TDim, TDim, vecs[ivec].begin ());
  T8_ASSERT (vecs.size () == num_vectors);
  return vecs;
}

/** Test the vector/point views */
TEST (t8_gtest_vec, vec_view)
{
  constexpr size_t seed = 12345;
  constexpr double min = -1e10, max = 1e10;
  std::mt19937_64 rng (seed);
  std::uniform_real_distribution<double> dist (min, max);
  constexpr size_t num_points = 10;

  /* Test for each dimension */
  unrolled_for (1, 4, idim, {
    double c_vectors[idim * num_points] = { 0 };
    /* Fill test vectors and create views. */
    for (size_t icoord = 0; icoord < num_points * idim; ++icoord)
      c_vectors[icoord] = dist (rng);
    auto vecs = t8_convert_array_to_vec<idim> (c_vectors, num_points);
    auto vec_views = t8_convert_array_to_vec_view<idim> (c_vectors, num_points);

    for (size_t ipoint = 0; ipoint < num_points; ++ipoint) {
      EXPECT_VEC_EQ (vecs[ipoint], vec_views[ipoint], T8_PRECISION_SQRT_EPS);
      /* Also check if functions return the same, but only for 3D. */
      if constexpr (idim == 3) {
        /* Normalize c vectors and cpp vectors. */
        t8_normalize (c_vectors + ipoint * idim);
        t8_normalize (vecs[ipoint]);
        /* Copied vector and vector view should be the same. */
        EXPECT_VEC_EQ (vecs[ipoint], vec_views[ipoint], T8_PRECISION_SQRT_EPS);
      }
    }
  });
}
