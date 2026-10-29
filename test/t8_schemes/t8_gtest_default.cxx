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
#include <t8_schemes/t8_scheme.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <utility>
#include <variant>

TEST (gtest_default_scheme, is_default)
{
  /* TODO: Implement an EXPECT_FALSE check for non-default schemes */
  const t8_scheme *scheme = t8_scheme_new_default ();
  for (int eclass = T8_ECLASS_VERTEX; eclass < T8_ECLASS_COUNT; ++eclass) {
    EXPECT_TRUE (t8_eclass_scheme_is_default (scheme, static_cast<t8_eclass_t> (eclass)));
  }
  scheme->unref ();
}

/** The eclass scheme at position \a TIndex of \ref t8_scheme::scheme_var. */
template <size_t TIndex>
using t8_test_eclass_scheme = std::variant_alternative_t<TIndex, t8_scheme::scheme_var>;

/** Get the size of the element type of the eclass scheme that \a scheme uses for \a tree_class.
 * \param [in] scheme     The scheme.
 * \param [in] tree_class The eclass of the current tree.
 * \tparam TIndices       The indices of all variant alternatives.
 * \return                The size of the element type, or 0 if no alternative of
 *                        \ref t8_scheme::scheme_var matches.
 */
template <size_t... TIndices>
static size_t
t8_test_get_element_type_size (const t8_scheme *scheme, const t8_eclass_t tree_class, std::index_sequence<TIndices...>)
{
  /** The size of the element type, 0 as long as no matching eclass scheme was found. */
  size_t element_type_size = 0;
  ((scheme->check_eclass_scheme_type<t8_test_eclass_scheme<TIndices>> (tree_class)
    && (element_type_size = sizeof (typename t8_test_eclass_scheme<TIndices>::element_type), true))
   || ...);
  return element_type_size;
}

struct element_type: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (GetParam ());
  }
  void
  TearDown () override
  {
    scheme->unref ();
  }

  const t8_scheme *scheme;
  t8_eclass_t eclass;
};

/* The element size of each eclass scheme must match the size of its element type. */
TEST_P (element_type, size_matches_element_size)
{
  const size_t element_type_size = t8_test_get_element_type_size (
    scheme, eclass, std::make_index_sequence<std::variant_size_v<t8_scheme::scheme_var>> {});
  EXPECT_NE (element_type_size, 0u);
  EXPECT_EQ (scheme->get_element_size (eclass), element_type_size);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_type, element_type, AllSchemes, print_all_schemes);
