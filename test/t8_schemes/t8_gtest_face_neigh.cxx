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
#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_element_c_interface.h>
#include <test/t8_gtest_macros.hxx>

#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>

/* *INDENT-OFF* */
class face_neigh: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    scheme = t8_scheme_new_default_cxx ();

    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &neigh);
    ts->t8_element_root (element);
  }

  void
  TearDown () override
  {
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &neigh);
    t8_scheme_cxx_unref (&scheme);
  }
  t8_element_t *element;
  t8_element_t *child;
  t8_element_t *neigh;
  t8_scheme_cxx *scheme;
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;

#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif
};
/* *INDENT-ON* */

void
t8_test_face_neighbor_inside (int num_faces, t8_element_t *element, t8_element_t *child, t8_element_t *neigh,
                              t8_eclass_scheme_c *ts)
{
  int face_num;
  int check;

  for (int iface = 0; iface < num_faces; iface++) {
    /* Compute the neighbors neighbor along a given face and check, if the result is the
     * original element. */
    ts->t8_element_face_neighbor_inside (child, neigh, iface, &face_num);
    ts->t8_element_face_neighbor_inside (neigh, element, face_num, &check);

    EXPECT_TRUE (ts->t8_element_equal (child, element)) << "Got a false neighbor.";
    EXPECT_ELEM_EQ (ts, child, element);
  }
}

int
t8_test_get_middle_child (t8_eclass_t eclass, int ilevel, t8_element_t *element, t8_element_t *child,
                          t8_eclass_scheme_c *ts)
{
  /* Get the child number of the child in the middle of the element, depending of the shape of the element. */
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return 0;
  case T8_ECLASS_LINE:
    return 0;
  case T8_ECLASS_QUAD:
    /* There are no inner children in level one refinement. The test starts with level two, because this is the first level, inner children exists.
       The third child of level one child 0 is one of four middle children in level two. */
    ts->t8_element_child (element, 0, child);
    ts->t8_element_copy (child, element);
    return 3;
  case T8_ECLASS_TRIANGLE:
    return 3;
  case T8_ECLASS_HEX:
    /* There are no inner children in level one refinement. The test starts with level two, because this is the first level, inner children existing.
       The third child of level one child 4 is one of eight middle children in level two. */
    ts->t8_element_child (element, 4, child);
    ts->t8_element_copy (child, element);
    return 3;
  case T8_ECLASS_TET:
    return 3;
  case T8_ECLASS_PRISM:
    /* There are no inner children in level one refinement. The test starts with level two, because this is the first level, inner children existing.
       The last child of level one child 4 is one of eight middle children in level two. */
    ts->t8_element_child (element, 4, child);
    ts->t8_element_copy (child, element);
    return 7;
  case T8_ECLASS_PYRAMID: {
    t8_dpyramid_t *pyramid = (t8_dpyramid_t *) element;
    /* middle_child_id of Pyramid Type 6. */
    if (pyramid->pyramid.type == T8_DPYRAMID_FIRST_TYPE) {
      return 8;
    }
    /* middle_child_id of Pyramid Type 7. */
    else {
      return 3;
    }
  }
  default:
    return 0;
  }
}

/* Compute all children along all faces. Compute their neighbors along the face,
 * check, if the children have root contact, and if the neighbors are outside of the
 * root. */
TEST_P (face_neigh, check_not_inside_root)
{
  /* Are the neighbors of the element really outside?. */

  const int num_faces = ts->t8_element_num_faces (element);

  for (int iface = 0; iface < num_faces; iface++) {

    const int num_children = ts->t8_element_num_face_children (element, iface);
    int *child_indices = T8_ALLOC (int, num_children);
    t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
    ts->t8_element_new (num_children, children);
    ts->t8_element_children_at_face (element, iface, children, num_children, child_indices);

    for (int jchild = 0; jchild < num_children; jchild++) {

      const int child_id = child_indices[jchild];
      const int face_contact = ts->t8_element_face_child_face (element, iface, jchild);

      ts->t8_element_child (element, child_id, child);
      int face_num;
      int inside = ts->t8_element_face_neighbor_inside (child, neigh, face_contact, &face_num);

      ASSERT_EQ (inside, 0) << "Element is not outside.";

      inside = ts->t8_element_tree_face (child, face_contact);
      ASSERT_EQ (inside, iface) << "Wrong face.";
    }
    ts->t8_element_destroy (num_children, children);
    T8_FREE (children);
    T8_FREE (child_indices);
  }
}

void
t8_recursive_check_diff (t8_element_t *element, t8_element_t *child, t8_element_t *neigh, t8_eclass_scheme_c *ts,
                         int maxlvl, int level)
{

  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  if (level == maxlvl) {
    return;
  }

  /* Compute the neighbors neighbor along a given face and check, if the result is the
   * original element. */
  int num_faces = ts->t8_element_num_faces (element);

  t8_test_face_neighbor_inside (num_faces, child, element, neigh, ts);

  int num_children = ts->t8_element_num_children (child);
  for (int ichild = 0; ichild < num_children; ichild++) {
    ts->t8_element_child (element, ichild, child);
    t8_recursive_check_diff (child, element, neigh, ts, maxlvl, level + 1);
    ts->t8_element_parent (child, element);
  }
}

/* Recursively check, if all neighbors are computed correct up to a given level. */
TEST_P (face_neigh, recursive_check_diff)
{
  int level = 1;
  int middle_child_id = t8_test_get_middle_child (eclass, level, element, child, ts);
  ts->t8_element_child (element, middle_child_id, child);

  t8_recursive_check_diff (child, element, neigh, ts, maxlvl, level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neigh, face_neigh, AllEclasses, print_eclass);
