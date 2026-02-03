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

/** \file t8_gtest_custom_assertion.hxx
* Provide customized GoogleTest functions for improved error-output
*/

#ifndef T8_GTEST_CUSTOM_ASSERTION_HXX
#define T8_GTEST_CUSTOM_ASSERTION_HXX

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>
#include <t8_forest/t8_forest_general.h>

/**
 * \brief Test two elements for equality and print the elements if they aren't equal
 *
 * \param[in] ts_expr The name of the scheme \a ts
 * \param[in] tree_class_expr The name of the tree class
 * \param[in] elem_1_expr The name of the first element \a elem_1
 * \param[in] elem_2_expr The name of the second element \a elem_2
 * \param[in] scheme The scheme to use to check the equality
 * \param[in] tree_class The eclass of the tree the elements are part of
 * \param[in] elem_1 The element to compare with \a elem_2
 * \param[in] elem_2 the element to compare with \a elem_1
 * \return testing::AssertionResult
 */
testing::AssertionResult
element_equality (const char *ts_expr, const char *tree_class_expr, const char *elem_1_expr, const char *elem_2_expr,
                  const t8_scheme *scheme, const t8_eclass_t eclass, const t8_element_t *elem_1,
                  const t8_element_t *elem_2)
{
  if (scheme->element_is_equal (eclass, elem_1, elem_2)) {
    return testing::AssertionSuccess (false);
  }
  else {
#if T8_ENABLE_DEBUG
    char elem_1_string[BUFSIZ];
    char elem_2_string[BUFSIZ];
    const t8_eclass_t tree_class = scheme->get_eclass_scheme_eclass (eclass);
    scheme->element_to_string (eclass, elem_1, elem_1_string, BUFSIZ);
    scheme->element_to_string (eclass, elem_2, elem_2_string, BUFSIZ);
    return testing::AssertionFailure (false)
           << elem_1_expr << " " << elem_1_string << " is not equal to \n"
           << elem_2_expr << " " << elem_2_string << " given scheme " << ts_expr << " and tree class "
           << tree_class_expr << " " << t8_eclass_to_string[tree_class];
#else
    return testing::AssertionFailure (false)
           << elem_1_expr << " is not equal to \n"
           << elem_2_expr << " given scheme " << ts_expr << " and tree class " << tree_class_expr;
#endif
  }
}

#define EXPECT_ELEM_EQ(scheme, eclass, elem1, elem2) \
  EXPECT_PRED_FORMAT4 (element_equality, (scheme), (eclass), (elem1), (elem2))

#define ASSERT_ELEM_EQ(scheme, tree_class, elem1, elem2) \
  ASSERT_PRED_FORMAT4 (element_equality, (scheme), (tree_class), (elem1), (elem2))

/**
 * \brief Test if two 3D Dimensionaltors are equal with respect to a given precision
 *
 * \tparam TDimensional2 Type of the first Dimensionaltor.
 * \param[in] Dimensional_1_expr Name of the first Dimensionaltor
 * \param[in] Dimensional_2_expr Name of the second Dimensionaltor
 * \param[in] precision_expr Name of the precision
 * \param[in] Dimensional_1 First Dimensionaltor to compare
 * \param[in] Dimensional_2 Second Dimensionaltor to compare
 * \param[in] precision Test equality up to this precision
 * \return testing::AssertionResult
 */
template <T8DimensionalType TDimensional1, T8DimensionalType TDimensional2>
static inline testing::AssertionResult
dimensional_equality (const char *Dimensional_1_expr, const char *Dimensional_2_expr, const char *precision_expr,
                      const TDimensional1 &Dimensional_1, const TDimensional2 &Dimensional_2, const double precision)
{
  if (t8_eq (Dimensional_1, Dimensional_2, precision)) {
    return testing::AssertionSuccess ();
  }
  else {
    return testing::AssertionFailure () << Dimensional_1_expr << " is not equal to " << Dimensional_2_expr << " \n"
                                        << "Precision given by " << precision_expr << " " << precision;
  }
}

#define EXPECT_VEC_EQ(Dimensional_1, Dimensional_2, precision) \
  EXPECT_PRED_FORMAT3 (dimensional_equality, (Dimensional_1), (Dimensional_2), (precision))

/**
 * \brief Test two forests for equality.
 *
 * \param[in] forest_A_expr The name of the forest \a forest_A
 * \param[in] forest_B_expr The name of the forest \a forest_B
 * \param[in] forest_A      The forest to compare with \a forest_B
 * \param[in] forest_B      The forest to compare with \a forest_A
 * \return testing::AssertionResult
 */
testing::AssertionResult
forest_equality (const char *forest_A_expr, const char *forest_B_expr, const t8_forest_t forest_A,
                 const t8_forest_t forest_B)
{
  if (t8_forest_is_equal (forest_A, forest_B)) {
    return testing::AssertionSuccess ();
  }
  else {
    return testing::AssertionFailure () << forest_A_expr << " is not equal to " << forest_B_expr;
  }
}

#define EXPECT_FOREST_EQ(forest_A, forest_B) EXPECT_PRED_FORMAT2 (forest_equality, (forest_A), (forest_B))

#endif /* T8_GTEST_CUSTOM_ASSERTION_HXX */
