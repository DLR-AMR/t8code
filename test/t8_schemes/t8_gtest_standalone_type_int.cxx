/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8.h"
#include "t8_element.h"
#include "t8_gtest_dfs_base.hxx"
#include "t8_schemes/t8_standalone/t8_standalone_elements.hxx"
#include "t8_schemes/t8_standalone/t8_standalone_implementation.hxx"


/* t8_standalone_scheme declares this class a friend, so that it can access determine_type_int. */
struct class_test_equal: public TestDFS
{
 private:
  /** Compute the type int of an element for the type bit \a itype manually: the bit at position maxlevel - level is
   * the type bit of the ancestor of the element at that level. For the levels finer than the element it is the type bit
   * of the element itself, which is the type bit of its first descendant at that level. */
  template <t8_eclass_t TEclass>
  static t8_element_coord
  manually_compute_type_int (const t8_standalone<TEclass> *element, const int itype)
  {
    using scheme = t8_standalone_scheme<TEclass>;
    const int maxlevel = scheme::get_maxlevel ();
    t8_standalone<TEclass> anc;
    scheme::element_copy ((const t8_element_t *) element, (t8_element_t *) &anc);
    t8_element_coord type_int = 0;
    for (int level = maxlevel; level >= anc.level; level--) {
      type_int |= (t8_element_coord) scheme::get_typebit (anc.type, itype) << (maxlevel - level);
    }
    while (anc.level) {
      scheme::element_get_parent ((const t8_element_t *) &anc, (t8_element_t *) &anc);
      type_int |= (t8_element_coord) scheme::get_typebit (anc.type, itype) << (maxlevel - anc.level);
    }
    return type_int;
  }

  template <t8_eclass_t TEclass>
  static void
  check_type_ints (const t8_element_t *element)
  {
    t8_standalone<TEclass> el;
    t8_standalone_scheme<TEclass>::element_copy (element, (t8_element_t *) &el);
    for (int itype = 0; itype < T8_ELEMENT_NUM_EQUATIONS[TEclass]; itype++) {
      EXPECT_EQ (manually_compute_type_int<TEclass> (&el, itype),
                 (t8_element_coord) t8_standalone_scheme<TEclass>::determine_type_int (&el, itype))
        << "itype " << itype << " level " << (int) el.level;
    }
  }

  void
  check_element () override
  {
    switch (eclass) {
    case T8_ECLASS_TRIANGLE:
      check_type_ints<T8_ECLASS_TRIANGLE> (element);
      break;
    case T8_ECLASS_TET:
      check_type_ints<T8_ECLASS_TET> (element);
      break;
    case T8_ECLASS_PRISM:
      check_type_ints<T8_ECLASS_PRISM> (element);
      break;
    case T8_ECLASS_PYRAMID:
      check_type_ints<T8_ECLASS_PYRAMID> (element);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }
};

TEST_P (class_test_equal, test_equal_dfs)
{
  check_recursive_dfs_to_max_lvl (5);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_standalone_type_int, class_test_equal,
                          ::testing::Combine (::testing::Values (1),
                                              ::testing::Values (T8_ECLASS_TRIANGLE, T8_ECLASS_TET, T8_ECLASS_PRISM,
                                                                 T8_ECLASS_PYRAMID)),
                          print_all_schemes);
