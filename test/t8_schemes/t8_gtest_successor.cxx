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
#include <t8_gtest_schemes.hxx>
#include <t8_gtest_custom_assertion.hxx>
#include <t8_gtest_macros.hxx>

class class_successor: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (GetParam ());
    scheme = create_from_scheme_id (scheme_id);
    tree_class = std::get<1> (GetParam ());
    scheme->element_new (tree_class, 1, &element);
    scheme->element_new (tree_class, 1, &successor);
    scheme->element_new (tree_class, 1, &child);
    scheme->element_new (tree_class, 1, &last);

    scheme->set_to_root (tree_class, element);

    tree_class = scheme->get_eclass_scheme_eclass (tree_class);
    if (tree_class == T8_ECLASS_VERTEX)
      GTEST_SKIP ();
  }
  void
  TearDown () override
  {
    scheme->element_destroy (tree_class, 1, &element);
    scheme->element_destroy (tree_class, 1, &successor);
    scheme->element_destroy (tree_class, 1, &child);
    scheme->element_destroy (tree_class, 1, &last);
    scheme->unref ();
  }
  const t8_scheme *scheme;
  t8_element_t *element;
  t8_element_t *successor;
  t8_element_t *child;
  t8_element_t *last;
  t8_eclass_t tree_class;
};

/* Check the computation of the successor recursively. Iterate through the elements
 * via DFS. On the given maxlvl-1 the children are computeted iteratively. For
 * each child, the successor is checked.
 */
static void
t8_recursive_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, t8_element_t *last,
                        const t8_scheme *scheme, const t8_eclass_t tree_class, const int maxlvl)
{
  const int level = scheme->element_get_level (tree_class, element);
  ASSERT_TRUE (scheme->element_get_level (tree_class, element) <= maxlvl
               && maxlvl <= scheme->get_maxlevel (tree_class) - 1);
  const int num_children = scheme->element_get_num_children (tree_class, element);
  if (level == maxlvl - 1) {
    /* Check, if the successor of the last recursion is the first child of
     * of this element.
     */
    scheme->element_get_child (tree_class, element, 0, child);
    EXPECT_ELEM_EQ (scheme, tree_class, child, successor);
    /*Check if the successor in this element is computed correctly */
    for (int ichild = 1; ichild < num_children; ichild++) {
      EXPECT_EQ (scheme->element_get_level (tree_class, child), maxlvl);
      scheme->element_construct_successor (tree_class, child, successor);
      scheme->element_get_child (tree_class, element, ichild, child);
      EXPECT_ELEM_EQ (scheme, tree_class, child, successor);
    }
    /*If the iterator is the last element, the test can finish */
    if (scheme->element_is_equal (tree_class, last, child)) {
      return;
    }
    /*Compute the next successor / "jump" out of the current element */
    else {
      EXPECT_EQ (scheme->element_get_level (tree_class, child), maxlvl);
      scheme->element_construct_successor (tree_class, child, successor);
    }
  }
  else {
    /*DFS run through the elements */
    for (int ichild = 0; ichild < num_children; ichild++) {
      scheme->element_get_child (tree_class, element, ichild, child);
      t8_recursive_successor (child, successor, element, last, scheme, tree_class, maxlvl);
      scheme->element_get_parent (tree_class, child, element);
    }
  }
}

/* Check the computation of the successor at the maximum level of the element.
 * Given the element of level maxlevel-2  with linear id 0 all children of
 * maximum level are computed. The successor runs through all these children.
 */
static void
t8_deep_successor (t8_element_t *element, t8_element_t *successor, t8_element_t *child, const t8_scheme *scheme,
                   const t8_eclass_t tree_class)
{
  const int maxlvl = scheme->get_maxlevel (tree_class);
  const int num_children = scheme->element_get_num_children (tree_class, element);

  for (int ichild = 0; ichild < num_children; ichild++) {
    scheme->element_get_child (tree_class, element, ichild, child);
    /* Go to the children at maximum level. */
    const int num_children_child = scheme->element_get_num_children (tree_class, child);
    for (int jchild = 0; jchild < num_children_child; jchild++) {
      scheme->element_get_child (tree_class, child, jchild, element);
      /* Check the computation of the successor. */
      ASSERT_ELEM_EQ (scheme, tree_class, element, successor) << "Wrong Successor at Maxlvl.\n";
      /* Compute the next successor. */
      EXPECT_EQ (scheme->element_get_level (tree_class, successor), maxlvl);
      scheme->element_construct_successor (tree_class, successor, successor);
    }
    scheme->element_get_parent (tree_class, child, element);
  }
}

TEST_P (class_successor, test_recursive_and_deep_successor)
{
#if T8_TEST_LEVEL_INT >= 1
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif

  /* Test at lower level. */
  for (int ilevel = 1; ilevel <= maxlvl; ilevel++) {
    scheme->element_set_linear_id (tree_class, successor, ilevel, 0);
    scheme->element_get_last_descendant (tree_class, element, last, ilevel);
    t8_recursive_successor (element, successor, child, last, scheme, tree_class, ilevel);
  }
  /* Test at Maxlevel. */
  scheme->element_set_linear_id (tree_class, element, scheme->get_maxlevel (tree_class) - 2, 0);
  scheme->element_set_linear_id (tree_class, successor, scheme->get_maxlevel (tree_class), 0);
  t8_deep_successor (element, successor, last, scheme, tree_class);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_successor, class_successor, AllSchemes, print_all_schemes);
