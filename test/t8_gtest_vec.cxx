/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/* In this file we collect tests for the routines in t8_vec.h */

#include <gtest/gtest.h>
#include <t8_vec.h>

/* Wrapper for 3D vector dataype */
typedef double t8_test_vec[3];
/* Accuracy used for comparisons with correct result */
#define epsilon 1e-9

/* test the t8_vec_norm function */
TEST (t8_gtest_vec, norm) {
  const t8_test_vec zero = {0, 0, 0};
  const t8_test_vec onetwothree = {1, 2, 3};
  const t8_test_vec arbitrary = {-.05, 3.14159, 42};

  const double normonetwothree = sqrt (1 + 4 + 9);
  const double normarbitrary = 42.117360883;

  EXPECT_EQ (t8_vec_norm (zero), 0);
  EXPECT_NEAR (t8_vec_norm (onetwothree), normonetwothree, epsilon);
  EXPECT_NEAR (t8_vec_norm (arbitrary), normarbitrary, epsilon);
}

/* test the t8_vec_dist function */
TEST (t8_gtest_vec, dist) {
  const t8_test_vec zero = {0, 0, 0};
  const t8_test_vec onetwothree = {1, 2, 3};
  const t8_test_vec arbitrary = {-.05, 3.14159, 42};
  const double distzeroonetwothree = sqrt (1 + 4 + 9);
  const double distarbitraryonetwothree = 39.030830477;

  EXPECT_EQ (t8_vec_dist (zero, zero), 0);
  EXPECT_EQ (t8_vec_dist (onetwothree, onetwothree), 0);
  EXPECT_NEAR (t8_vec_dist (onetwothree, zero), distzeroonetwothree, epsilon);
  EXPECT_NEAR (t8_vec_dist (zero, onetwothree), distzeroonetwothree, epsilon);
  EXPECT_NEAR (t8_vec_dist (arbitrary, onetwothree), distarbitraryonetwothree, epsilon);
  EXPECT_NEAR (t8_vec_dist (onetwothree, arbitrary), distarbitraryonetwothree, epsilon);
}

/* test the t8_vec_ax function */
TEST (t8_gtest_vec, ax) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec conetwothree = {1, 2, 3};
  const t8_test_vec carbitrary = {-.05, 3.14159, 42};
  t8_test_vec zero = {0, 0, 0};
  t8_test_vec onetwothree = {1, 2, 3};
  t8_test_vec arbitrary = {-.05, 3.14159, 42};
  const double alpha = 5.1234;

  /* Compute Y = alpha * Y */
  t8_vec_ax (zero, alpha);
  t8_vec_ax (onetwothree, alpha);
  t8_vec_ax (arbitrary, alpha);

  /* Check results */
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (zero[i], alpha * czero[i], epsilon);
    EXPECT_NEAR (onetwothree[i], alpha * conetwothree[i], epsilon);
    EXPECT_NEAR (arbitrary[i], alpha * carbitrary[i], epsilon);
  }
}

/* test the t8_vec_axb function */
TEST (t8_gtest_vec, axb) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec conetwothree = {1, 2, 3};
  const t8_test_vec carbitrary = {-.05, 3.14159, 42};
  const double b = 2.71828;
  t8_test_vec zero;
  t8_test_vec onetwothree;
  t8_test_vec arbitrary;
  const double alpha = 5.1234;

  /* Compute Y = alpha * Y + b */
  t8_vec_axb (czero, zero, alpha, b);
  t8_vec_axb (conetwothree, onetwothree, alpha, b);
  t8_vec_axb (carbitrary, arbitrary, alpha, b);

  /* Check results */
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (zero[i], alpha * czero[i] + b, epsilon);
    EXPECT_NEAR (onetwothree[i], alpha * conetwothree[i] + b, epsilon);
    EXPECT_NEAR (arbitrary[i], alpha * carbitrary[i] + b, epsilon);
  }
}

/* test the t8_vec_axpy function */
TEST (t8_gtest_vec, axpy) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec conetwothree = {1, 2, 3};
  const t8_test_vec carbitrary = {-.05, 3.14159, 42};
  const t8_test_vec init = {3, 2.71828, -4.1};
  /* The next three vecs must be initialized with the values of init */
  t8_test_vec zero = {3, 2.71828, -4.1};
  t8_test_vec onetwothree = {3, 2.71828, -4.1};
  t8_test_vec arbitrary {3, 2.71828, -4.1};
  const double alpha = 5.1234;

  /* Compute Y = Y + alpha * Y */
  t8_vec_axpy (czero, zero, alpha);
  t8_vec_axpy (conetwothree, onetwothree, alpha);
  t8_vec_axpy (carbitrary, arbitrary, alpha);

  /* Check results */
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (zero[i], init[i] + alpha * czero[i], epsilon);
    EXPECT_NEAR (onetwothree[i], init[i] + alpha * conetwothree[i], epsilon);
    EXPECT_NEAR (arbitrary[i], init[i] + alpha * carbitrary[i], epsilon);
  }
}

/* test the t8_vec_axpyz function */
TEST (t8_gtest_vec, axpyz) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec conetwothree = {1, 2, 3};
  const t8_test_vec carbitrary = {-.05, 3.14159, 42};
  const t8_test_vec init = {3, 2.71828, -4.1};
  t8_test_vec Z;
  const double alpha = 5.1234;

  /* Z = init + alpha * zero */
  t8_vec_axpyz (czero, init, Z, alpha);
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * czero[i], epsilon);
  }
  /* Z = init + alpha * conetwothree */
  t8_vec_axpyz (conetwothree, init, Z, alpha);
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * conetwothree[i], epsilon);
  }
  /* Z = init + alpha * carbitrary */
  t8_vec_axpyz (carbitrary, init, Z, alpha);
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (Z[i], init[i] + alpha * carbitrary[i], epsilon);
  }
}

/* test the t8_vec_dot function */
TEST (t8_gtest_vec, dot) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec conetwothree = {1, 2, 3};
  const t8_test_vec carbitrary = {-.05, 3.14159, 42};
  double result;

  /* Dot product with 0 is 0 */
  for (int i = 0;i < 3;++i) {
    EXPECT_EQ (t8_vec_dot (czero, czero), 0);
    EXPECT_EQ (t8_vec_dot (czero, conetwothree), 0);
    EXPECT_EQ (t8_vec_dot (czero, carbitrary), 0);
  }

  /* Test dot product of the remaining two vecs */
  result = 0;
  for (int i = 0;i < 3;++i) {
    result += conetwothree[i] * carbitrary[i];
  }
  EXPECT_NEAR (t8_vec_dot (conetwothree, carbitrary), result, epsilon);

  /* For the dot-product of a vector with itself we use the square of its norm */
  result = t8_vec_norm (conetwothree) * t8_vec_norm (conetwothree);
  EXPECT_NEAR (t8_vec_dot (conetwothree, conetwothree), result, epsilon);

  result = t8_vec_norm (carbitrary) * t8_vec_norm (carbitrary);
  EXPECT_NEAR (t8_vec_dot (carbitrary, carbitrary), result, epsilon);
}

/* test the t8_vec_cross function */
TEST (t8_gtest_vec, cross) {
  const t8_test_vec czero = {0, 0, 0};
  const t8_test_vec e1 = {1, 0, 0};
  const t8_test_vec e2 = {0, 1, 0};
  const t8_test_vec e3 = {0, 0, 1};
  t8_test_vec cross;

  /* cross product with 0 is 0 */
  t8_vec_cross (czero, czero, cross);
  for (int i = 0;i < 3;++i) {
    EXPECT_EQ (cross[i], 0);
  }
  t8_vec_cross (e1, czero, cross);
  for (int i = 0;i < 3;++i) {
    EXPECT_EQ (cross[i], 0);
  }
  t8_vec_cross (e2, czero, cross);
  for (int i = 0;i < 3;++i) {
    EXPECT_EQ (cross[i], 0);
  }

  /* e1 x e2 = e3 */
  t8_vec_cross (e1, e2, cross);
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (cross[i], e3[i], epsilon);
  }
  /* e2 x e3 = e1 */
  t8_vec_cross (e2, e3, cross);
  for (int i = 0;i < 3;++i) {
    EXPECT_NEAR (cross[i], e1[i], epsilon);
  }
}
