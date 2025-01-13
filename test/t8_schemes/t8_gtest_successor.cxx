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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>

class class_successor: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    scheme = t8_scheme_all_schemes ();
    scheme_id = GetParam ();
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &element);
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &successor);
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &child);
    scheme->element_new (static_cast<t8_eclass_t> (scheme_id), 1, &last);

    scheme->get_root (static_cast<t8_eclass_t> (scheme_id), element);

    tree_class = scheme->get_eclass_scheme_eclass (static_cast<t8_eclass_t> (scheme_id));
    if (tree_class == T8_ECLASS_VERTEX)
      GTEST_SKIP ();
  }
  void
  TearDown () override
  {
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &element);
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &successor);
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &child);
    scheme->element_destroy (static_cast<t8_eclass_t> (scheme_id), 1, &last);
    scheme->unref ();
  }
  t8_eclass_t tree_class;
  t8_scheme *scheme;
  t8_element_t *element;
  t8_element_t *successor;
  t8_element_t *child;
  t8_element_t *last;
};

/* Check the computation of the successor recursively. Iterate through the elements
 * via DFS. On the given maxlvl-1 the children are computeted iteratively. For
 * each child, the successor is checked.
 */
static void
t8_recursive_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, t8_element_t *last,
                        t8_scheme *scheme, const int scheme_id, const int maxlvl)
{
  const int level = scheme->element_get_level (static_cast<t8_eclass_t> (scheme_id), element);
  ASSERT_TRUE (scheme->element_get_level (static_cast<t8_eclass_t> (scheme_id), element) <= maxlvl
               && maxlvl <= scheme->get_maxlevel (static_cast<t8_eclass_t> (scheme_id)) - 1);
  const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t> (scheme_id), element);
  if (level == maxlvl - 1) {
    /* Check, if the successor of the last recursion is the first child of
     * of this element.
     */
    scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, 0, child);
    EXPECT_ELEM_EQ (scheme, scheme_id, child, successor);
    /*Check if the successor in this element is computed correctly */
    for (int ichild = 1; ichild < num_children; ichild++) {
      EXPECT_EQ (scheme->element_get_level (static_cast<t8_eclass_t> (scheme_id), child), maxlvl);
      scheme->element_construct_successor (static_cast<t8_eclass_t> (scheme_id), child, successor);
      scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, ichild, child);
      EXPECT_ELEM_EQ (scheme, scheme_id, child, successor);
    }
    /*If the iterator is the last element, the test can finish */
    if (scheme->element_is_equal (static_cast<t8_eclass_t> (scheme_id), last, child)) {
      return;
    }
    /*Compute the next successor / "jump" out of the current element */
    else {
      EXPECT_EQ (scheme->element_get_level (static_cast<t8_eclass_t> (scheme_id), child), maxlvl);
      scheme->element_construct_successor (static_cast<t8_eclass_t> (scheme_id), child, successor);
    }
  }
  else {
    /*DFS run through the elements */
    for (int ichild = 0; ichild < num_children; ichild++) {
      scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, ichild, child);
      t8_recursive_successor (child, successor, element, last, scheme, scheme_id, maxlvl);
      scheme->element_get_parent (static_cast<t8_eclass_t> (scheme_id), child, element);
    }
  }
}

/* Check the computation of the successor at the maximum level of the element.
 * Given the element of level maxlevel-2  with linear id 0 all children of
 * maximum level are computed. The successor runs through all these children.
 */
static void
t8_deep_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, t8_scheme *scheme,
                   const int scheme_id)
{
  const int maxlvl = scheme->get_maxlevel (static_cast<t8_eclass_t> (scheme_id));
  const int num_children = scheme->element_get_num_children (static_cast<t8_eclass_t> (scheme_id), element);

  for (int ichild = 0; ichild < num_children; ichild++) {
    scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), element, ichild, child);
    /* Go to the children at maximum level. */
    const int num_children_child = scheme->element_get_num_children (static_cast<t8_eclass_t> (scheme_id), child);
    for (int jchild = 0; jchild < num_children_child; jchild++) {
      scheme->element_get_child (static_cast<t8_eclass_t> (scheme_id), child, jchild, element);
      /* Check the computation of the successor. */
      ASSERT_TRUE (scheme->element_is_equal (static_cast<t8_eclass_t> (scheme_id), element, successor))
        << "Wrong Successor at Maxlvl.\n";
      /* Compute the next successor. */
      EXPECT_EQ (scheme->element_get_level (static_cast<t8_eclass_t> (scheme_id), successor), maxlvl);
      scheme->element_construct_successor (static_cast<t8_eclass_t> (scheme_id), successor, successor);
    }
    scheme->element_get_parent (static_cast<t8_eclass_t> (scheme_id), child, element);
  }
}

TEST_P (class_successor, test_recursive_and_deep_successor)
{
#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  /* Test at lower level. */
  for (int ilevel = 1; ilevel <= maxlvl; ilevel++) {
    scheme->element_set_linear_id (static_cast<t8_eclass_t> (scheme_id), successor, ilevel, 0);
    scheme->element_get_last_descendant (static_cast<t8_eclass_t> (scheme_id), element, last, ilevel);
    t8_recursive_successor (element, successor, child, last, scheme, scheme_id, ilevel);
  }
  /* Test at Maxlevel. */
  scheme->element_set_linear_id (static_cast<t8_eclass_t> (scheme_id), element,
                                 scheme->get_maxlevel (static_cast<t8_eclass_t> (scheme_id)) - 2, 0);
  scheme->element_set_linear_id (static_cast<t8_eclass_t> (scheme_id), successor,
                                 scheme->get_maxlevel (static_cast<t8_eclass_t> (scheme_id)), 0);
  t8_deep_successor (element, successor, last, scheme, scheme_id);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_successor, class_successor, AllEclasses, print_eclass);
