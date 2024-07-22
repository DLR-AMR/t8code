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

/** \file t8_gtest_root.cxx
*/

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>

class root: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    scheme = t8_scheme_new_default_cxx ();
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &element);
    ts->t8_element_root (element);
  }
  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &element);
    t8_scheme_cxx_unref (&scheme);
  }
  t8_element_t *element;
  t8_scheme_cxx *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;
};

/*Test root*/

TEST_P (root, has_level_zero)
{
  EXPECT_EQ (ts->t8_element_level (element), 0);
}

TEST_P (root, equals_linear_id_0_0)
{
  t8_element_t *root_compare;
  ts->t8_element_new (1, &root_compare);
  ts->t8_element_set_linear_id (root_compare, 0, 0);
  EXPECT_ELEM_EQ (ts, element, root_compare);
  ts->t8_element_destroy (1, &root_compare);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_root, root, AllEclasses, print_eclass);
