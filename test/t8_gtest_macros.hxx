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
#include <test/t8_schemes/t8_gtest_iterator.cxx>
#include <iostream>

extern t8_scheme_cxx *
t8_scheme_new_default_cxx ();

/**
 * lambda to pass to an INSTANTIATE_TEST_SUITE_P to print the current cmesh_example_base
 * 
 */
auto print_eclass = [] (const testing::TestParamInfo<t8_eclass> &info) { return t8_eclass_to_string[info.param]; };

/**
 * Number of points to use in tests
 * 
 */
#ifdef T8_ENABLE_LESS_TESTS
#define T8_NUM_SAMPLE_POINTS 1000
#else
#define T8_NUM_SAMPLE_POINTS 10000
#endif

const t8_scheme_cxx *default_scheme = t8_scheme_new_default_cxx ();
const t8_scheme_cxx *sa_scheme = t8_scheme_new_default_cxx ();
const std::vector<const t8_scheme_cxx *> schemes = { default_scheme, sa_scheme };
scheme_iterators scheme_iter (schemes);

#define AllSchemesEclasses testing::ValuesIn (scheme_iter.begin (), scheme_iter.end ())
#define AllEclasses testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT)
#define AllEclasses2D testing::Values (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE)

#endif /* T8_GTEST_MACROS_HXX */
