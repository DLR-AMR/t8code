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

/** \file t8_gtest_default.cxx
* Provide tests to check the functionality of the default-scheme. 
*/

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_default_line.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid.hxx>
#include <test/t8_gtest_macros.hxx>

TEST (gtest_default_scheme, is_default)
{
  /* TODO: Implement an EXPECT_FALSE check for non-default schemes */
  t8_scheme *scheme = t8_scheme_new_default ();
  for (int eclass = T8_ECLASS_VERTEX; eclass < T8_ECLASS_COUNT; ++eclass) {
    EXPECT_TRUE (t8_eclass_scheme_is_default (scheme, static_cast<t8_eclass_t> (eclass)));
  }
  scheme->unref ();
}
