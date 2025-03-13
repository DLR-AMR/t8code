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

/** \file t8_gtest_macros.hxx
* Provide macros for instantiating parameterized tests
*/

#ifndef T8_GTEST_MACROS_HXX
#define T8_GTEST_MACROS_HXX

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <iostream>
#include <t8_schemes/t8_scheme.hxx>

/**
 * lambda to pass to an INSTANTIATE_TEST_SUITE_P to print the current cmesh_example_base
 * 
 */
inline auto print_eclass
  = [] (const testing::TestParamInfo<t8_eclass> &info) { return t8_eclass_to_string[info.param]; };

/**
 * Concept which checks if a class implements the function SetUpTestsuite ()
 * \tparam TTestingClass The class to check
 */
template <typename TTestingClass>
concept HasSetUpTestSuite = requires (TTestingClass testing_class) {
  {
    testing_class.SetUpTestSuite ()
  } -> std::same_as<void>;
};

/**
 * Testing class which registers a package so that the attribute system can be used.
 * \tparam TTestingClass The underlying parameterized testing class. Has to implement SetUpTestsuite ()
 */
template <class TTestingClass>
  requires HasSetUpTestSuite<TTestingClass>
class t8_test_with_attributes: public TTestingClass {
 public:
  /**
   * Registers a package id for the test system
   */
  static void
  SetUpTestSuite ()
  {
    TTestingClass::SetUpTestSuite ();
    t8_testsuite_package_id
      = sc_package_register (NULL, SC_LP_DEFAULT, ::testing::UnitTest::GetInstance ()->current_test_info ()->name (),
                             "t8code testsuite package. Used for testing of external user attributes.");
  }

  static int
  get_testsuite_package_id (void)
  {
    return t8_testsuite_package_id;
  }

 private:
  static int t8_testsuite_package_id; /** Package id of the testsuite */
};

/**
 * Number of points to use in tests
 * 
 */
#if T8CODE_TEST_LEVEL == 1
#define T8_NUM_SAMPLE_POINTS 1000
#else
#define T8_NUM_SAMPLE_POINTS 10000
#endif

#define AllEclasses testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT)
#define AllEclasses2D testing::Values (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE)

#endif /* T8_GTEST_MACROS_HXX */
