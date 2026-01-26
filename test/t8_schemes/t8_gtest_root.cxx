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
#include <test/t8_gtest_schemes.hxx>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

struct root: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
    scheme->element_new (eclass, 1, &element);
    scheme->set_to_root (eclass, element);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &element);
    scheme->unref ();
  }

  t8_element_t *element;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
};

/*Test root*/

TEST_P (root, has_level_zero)
{
  EXPECT_EQ (scheme->element_get_level (eclass, element), 0);
}

TEST_P (root, equals_linear_id_0_0)
{
  t8_element_t *root_compare;
  scheme->element_new (eclass, 1, &root_compare);
  scheme->element_set_linear_id (eclass, root_compare, 0, 0);
  EXPECT_ELEM_EQ (scheme, eclass, element, root_compare);
  scheme->element_destroy (eclass, 1, &root_compare);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_root, root, AllSchemes, print_all_schemes);
