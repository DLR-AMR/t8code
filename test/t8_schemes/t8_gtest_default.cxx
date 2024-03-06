/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 The University of Texas System
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

/** \file t8_gtest_default.cxx
* Provide tests to check the functionality of the default-scheme. 
*/

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

#include <t8_schemes/t8_default/t8_default_vertex/t8_default_vertex_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_line/t8_default_line_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_tet/t8_default_tet_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_prism/t8_default_prism_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid_cxx.hxx>
#include <test/t8_gtest_macros.hxx>

class gtest_default_scheme: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    /* Construct every eclass scheme explixitly */
    const t8_eclass_t eclass = GetParam ();
    switch (eclass) {
    case T8_ECLASS_VERTEX:
      eclass_scheme = new t8_default_scheme_vertex_c ();
      break;
    case T8_ECLASS_LINE:
      eclass_scheme = new t8_default_scheme_line_c ();
      break;
    case T8_ECLASS_QUAD:
      eclass_scheme = new t8_default_scheme_quad_c ();
      break;
    case T8_ECLASS_TRIANGLE:
      eclass_scheme = new t8_default_scheme_tri_c ();
      break;
    case T8_ECLASS_HEX:
      eclass_scheme = new t8_default_scheme_hex_c ();
      break;
    case T8_ECLASS_TET:
      eclass_scheme = new t8_default_scheme_tet_c ();
      break;
    case T8_ECLASS_PRISM:
      eclass_scheme = new t8_default_scheme_prism_c ();
      break;
    case T8_ECLASS_PYRAMID:
      eclass_scheme = new t8_default_scheme_pyramid_c ();
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }
  void
  TearDown () override
  {
    delete (eclass_scheme);
  }
  t8_eclass_scheme_c *eclass_scheme;
};

TEST_P (gtest_default_scheme, is_default)
{
  /* TODO: Implement an EXPECT_FALSE check for non-default schemes */
  EXPECT_TRUE (t8_eclass_scheme_is_default (eclass_scheme));
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_default_scheme, gtest_default_scheme, AllEclasses, print_eclass);
