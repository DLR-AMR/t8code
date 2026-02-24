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

/** \file t8_gtest_macros.hxx
* Provide macros for instantiating parameterized tests
*/

#ifndef T8_GTEST_MACROS_HXX
#define T8_GTEST_MACROS_HXX

#include <gtest/gtest.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <iostream>
#include <t8_schemes/t8_scheme.hxx>
#include <test/t8_gtest_memory_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

/**
 * lambda to pass to an INSTANTIATE_TEST_SUITE_P to print the current cmesh_example_base
 * 
 */
inline auto print_eclass
  = [] (const testing::TestParamInfo<t8_eclass> &info) { return t8_eclass_to_string[info.param]; };

/** Define a lambda to beautify gtest output for tuples <level, cmesh>.
 * This will set the correct level and cmesh name as part of the test case name. */
auto pretty_print_eclass_scheme_and_level
  = [] (const testing::TestParamInfo<std::tuple<std::tuple<int, t8_eclass_t>, int>> &info) {
      std::string scheme = t8_scheme_to_string[std::get<0> (std::get<0> (info.param))];
      std::string eclass = t8_eclass_to_string[std::get<1> (std::get<0> (info.param))];
      std::string level = std::string ("_level_") + std::to_string (std::get<1> (info.param));
      return scheme + "_" + eclass + level;
    };

/**
 * Initializes everything needed for the t8code testsuite.
 * MPI is initialized with MPI_COMM_WORLD and SC with loglevel SC_LP_PRODUCTION.
 * \param [in] argc Number of command-line arguments.
 * \param [in] argv Array of command-line argument strings. 
 * \param [in] log_threshold    The log threshold used to initialize t8code.
 */
void
t8_testsuite_init (int *argc, char ***argv, int log_threshold);

/**
 * Finalizes everything needed in the t8code testsuite.
 */
void
t8_testsuite_finalize ();

/**
 * Returns the attribute package id of the t8code testsuite.
 * \return The package id.
 */
int
t8_testsuite_get_package_id ();

/**
 * Number of points to use in tests
 * 
 */
#if T8_TEST_LEVEL_INT >= 2
#define T8_NUM_SAMPLE_POINTS 500
#elif T8_TEST_LEVEL_INT >= 1
#define T8_NUM_SAMPLE_POINTS 1000
#else
#define T8_NUM_SAMPLE_POINTS 10000
#endif

/** Google test range for parametrized covering all eclasses. */
#define AllEclasses testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT)
/** Google test range for parametrized covering all 2 dimensional eclasses. */
#define AllEclasses2D testing::Values (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE)

#endif /* T8_GTEST_MACROS_HXX */
