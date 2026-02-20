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

/** \file t8_gtest_schemes.hxx 
 * File to provide different schemes for testing purposes.
 */

#ifndef T8_GTEST_SCHEMES_HXX
#define T8_GTEST_SCHEMES_HXX

#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <gtest/gtest.h>

/** Create a scheme according to a scheme id.
 * \param [in] scheme_id 0: Use default scheme; 1: Use standalone scheme.
 * \return The created scheme.
 */
const t8_scheme *
create_from_scheme_id (const int scheme_id)
{
  switch (scheme_id) {
  case 0:
    return t8_scheme_new_default ();
  case 1:
    return t8_scheme_new_standalone ();
  default:
    SC_ABORT_NOT_REACHED ();
    return nullptr;
  }
}

/** Strings for the two scheme types. */
static const char *t8_scheme_to_string[] = { "default", "standalone" };

/** Lambda to print the scheme and the eclass of an TestParamInfo object. */
auto print_all_schemes = [] (const testing::TestParamInfo<std::tuple<int, t8_eclass_t>> &info) {
  return std::string (t8_scheme_to_string[std::get<0> (info.param)]) + "_"
         + t8_eclass_to_string[std::get<1> (info.param)];
};

/** Lambda to print the scheme. */
auto print_scheme
  = [] (const testing::TestParamInfo<int> &info) { return std::string (t8_scheme_to_string[info.param]); };

/** Macro for all schemes. */
#define AllSchemeCollections ::testing::Range (0, 2)
/** Macro for all schemes and all possible eclasses.*/
#define AllSchemes ::testing::Combine (AllSchemeCollections, ::testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT))

#endif /* T8_GTEST_SCHEMES_HXX */
