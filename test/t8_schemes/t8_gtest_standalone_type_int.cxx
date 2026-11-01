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
#include <vector>
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

  template <t8_eclass_t TEclass>
  static void
  collect_elements (const t8_standalone<TEclass> &element, const int maxlevel,
                    std::vector<t8_standalone<TEclass>> &elements)
  {
    using scheme = t8_standalone_scheme<TEclass>;
    elements.push_back (element);
    if (element.level < maxlevel) {
      const int num_children = scheme::element_get_num_children ((const t8_element_t *) &element);
      for (int ichild = 0; ichild < num_children; ichild++) {
        t8_standalone<TEclass> child;
        scheme::element_get_child ((const t8_element_t *) &element, ichild, (t8_element_t *) &child);
        collect_elements<TEclass> (child, maxlevel, elements);
      }
    }
  }

  /* Check that element_compare_type_int orders all pairs of elements up to \a maxlevel like element_compare. */
  template <t8_eclass_t TEclass>
  static void
  check_compare (const int maxlevel)
  {
    using scheme = t8_standalone_scheme<TEclass>;
    t8_standalone<TEclass> root;
    scheme::set_to_root ((t8_element_t *) &root);
    std::vector<t8_standalone<TEclass>> elements;
    collect_elements<TEclass> (root, maxlevel, elements);
    const auto sign = [] (const int value) { return (value > 0) - (value < 0); };
    size_t num_mismatches = 0;
    for (const auto &el1 : elements) {
      for (const auto &el2 : elements) {
        const int expected = sign (scheme::element_compare ((const t8_element_t *) &el1, (const t8_element_t *) &el2));
        const int result
          = sign (scheme::element_compare_type_int ((const t8_element_t *) &el1, (const t8_element_t *) &el2));
        if (expected != result && num_mismatches++ < 5) {
          ADD_FAILURE () << "levels " << (int) el1.level << ", " << (int) el2.level << " ids "
                         << scheme::element_get_linear_id ((const t8_element_t *) &el1, el1.level) << ", "
                         << scheme::element_get_linear_id ((const t8_element_t *) &el2, el2.level) << ": expected "
                         << expected << ", got " << result;
        }
      }
    }
    EXPECT_EQ (num_mismatches, 0u) << "of " << elements.size () * elements.size () << " pairs";
  }

  template <t8_eclass_t TEclass>
  static void
  check_element_of_class (const t8_element_t *element, const bool compare)
  {
    if (compare) {
      check_compare<TEclass> (T8_ELEMENT_DIM[TEclass] <= 2 ? 4 : 3);
    }
    else {
      check_type_ints<TEclass> (element);
    }
  }

  void
  check_class (const bool compare)
  {
    switch (eclass) {
    case T8_ECLASS_LINE:
      check_element_of_class<T8_ECLASS_LINE> (element, compare);
      break;
    case T8_ECLASS_QUAD:
      check_element_of_class<T8_ECLASS_QUAD> (element, compare);
      break;
    case T8_ECLASS_TRIANGLE:
      check_element_of_class<T8_ECLASS_TRIANGLE> (element, compare);
      break;
    case T8_ECLASS_HEX:
      check_element_of_class<T8_ECLASS_HEX> (element, compare);
      break;
    case T8_ECLASS_TET:
      check_element_of_class<T8_ECLASS_TET> (element, compare);
      break;
    case T8_ECLASS_PRISM:
      check_element_of_class<T8_ECLASS_PRISM> (element, compare);
      break;
    case T8_ECLASS_PYRAMID:
      check_element_of_class<T8_ECLASS_PYRAMID> (element, compare);
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }

  void
  check_element () override
  {
    check_class (false);
  }

 public:
  void
  check_compare_all_pairs ()
  {
    check_class (true);
  }
};

TEST_P (class_test_equal, test_equal_dfs)
{
  check_recursive_dfs_to_max_lvl (5);
}

TEST_P (class_test_equal, test_compare_type_int)
{
  check_compare_all_pairs ();
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_standalone_type_int, class_test_equal,
                          ::testing::Combine (::testing::Values (1),
                                              ::testing::Values (T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE,
                                                                 T8_ECLASS_HEX, T8_ECLASS_TET, T8_ECLASS_PRISM,
                                                                 T8_ECLASS_PYRAMID)),
                          print_all_schemes);
