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

/** \file t8_gtest_custom_assertion.cxx
* Provide customized GoogleTest functions for improved error-output
*/

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_vec.h>

#ifndef CUSTOM_ASSERTION_HXX
#define CUSTOM_ASSERTION_HXX

/**
 * \brief Test two elements for equality and print the elements if they aren't equal
 * 
 * \param[in] ts_expr The name of the scheme \a ts
 * \param[in] elem_1_expr The name of the first element \a elem_1
 * \param[in] elem_2_expr The name of the second element \a elem_2
 * \param[in] ts The scheme to use to check the equality
 * \param[in] elem_1 The element to compare with \a elem_2
 * \param[in] elem_2 the element to compare with \a elem_1
 * \return testing::AssertionResult 
 */
testing::AssertionResult
element_equality (const char *ts_expr, const char *elem_1_expr, const char *elem_2_expr, const t8_scheme_c *ts,
                  const t8_element_t *elem_1, const t8_element_t *elem_2)
{
  if (ts->t8_element_equal (elem_1, elem_2)) {
    return testing::AssertionSuccess ();
  }
  else {
#if T8_ENABLE_DEBUG
    char elem_1_string[BUFSIZ];
    char elem_2_string[BUFSIZ];
    ts->t8_element_to_string (elem_1, elem_1_string, BUFSIZ);
    ts->t8_element_to_string (elem_2, elem_2_string, BUFSIZ);
    return testing::AssertionFailure () << elem_1_expr << " " << elem_1_string << " is not equal to \n"
                                        << elem_2_expr << " " << elem_2_string << " given scheme " << ts_expr;
#else
    return testing::AssertionFailure () << elem_1_expr << " is not equal to \n"
                                        << elem_2_expr << " given scheme " << ts_expr;
#endif
  }
}

#define EXPECT_ELEM_EQ(scheme, elem1, elem2) EXPECT_PRED_FORMAT3 (element_equality, (scheme), (elem1), (elem2))

/**
 * \brief Test if two 3D vectors are equal with respect to a given precision
 * 
 * \param[in] vec_1_expr Name of the first vector
 * \param[in] vec_2_expr Name of the second vector
 * \param[in] precision_expr Name of the precision
 * \param[in] vec_1 First vector to compare
 * \param[in] vec_2 Second vector to compare
 * \param[in] precision Test equality up to this prescision
 * \return testing::AssertionResult 
 */
testing::AssertionResult
vec3_equality (const char *vec_1_expr, const char *vec_2_expr, const char *precision_expr, const double vec_1[3],
               const double vec_2[3], const double precision)
{
  if (t8_vec_eq (vec_1, vec_2, precision)) {
    return testing::AssertionSuccess ();
  }
  else {
#if T8_ENABLE_DEBUG
    return testing::AssertionFailure () << vec_1[0] << " " << vec_1[1] << " " << vec_1[2] << " " << vec_1_expr
                                        << " is not equal to \n"
                                        << vec_2[0] << " " << vec_2[1] << " " << vec_2[2] << " " << vec_2_expr << " \n"
                                        << "Precision given by " << precision_expr << " " << precision;
#else
    return testing::AssertionFailure () << vec_1_expr << " is not equal to " << vec_2_expr;
#endif
  }
}

#define EXPECT_VEC3_EQ(vec_1, vec_2, precision) EXPECT_PRED_FORMAT3 (vec3_equality, (vec_1), (vec_2), (precision))

#endif /* CUSTOM_ASSERTION_HXX */
