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

#ifndef T8_GTEST_SCHEMES_HXX
#define T8_GTEST_SCHEMES_HXX

#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_scheme_builder.hxx>
#include <gtest/gtest.h>

const t8_scheme *
create_from_scheme_id (const int scheme_id)
{
  switch (scheme_id) {
  case 0:
    return t8_scheme_new_default ();
  default:
    SC_ABORT_NOT_REACHED ();
    return nullptr;
  }
}

#define AllSchemes ::testing::Combine (::testing::Values (0), ::testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT))
#define AllSchemeCollections ::testing::Values (0)

#endif /* T8_GTEST_SCHEMES_HXX */
