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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_element_c_interface.h>
#include <test/t8_gtest_macros.hxx>

#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>

class face_neigh: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    tree_class = GetParam ();
    scheme = t8_scheme_new_default ();
    scheme->element_new (tree_class, 1, &element);
    scheme->element_new (tree_class, 1, &child);
    scheme->element_new (tree_class, 1, &neigh);
    scheme->get_root (tree_class, element);
  }

  void
  TearDown () override
  {
    scheme->element_destroy (tree_class, 1, &element);
    scheme->element_destroy (tree_class, 1, &child);
    scheme->element_destroy (tree_class, 1, &neigh);
    scheme->unref ();
  }
  t8_element_t *element;
  t8_element_t *child;
  t8_element_t *neigh;
  t8_scheme *scheme;
  t8_eclass_t tree_class;

#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 4;
#endif
};

void
t8_test_face_neighbor_inside (int num_faces, t8_element_t *element, t8_element_t *child, t8_element_t *neigh,
                              t8_scheme *scheme, t8_eclass_t tree_class)
{
  int face_num;
  int check;

  for (int iface = 0; iface < num_faces; iface++) {
    /* Compute the neighbors neighbor along a given face and check, if the result is the
     * original element. */
    scheme->element_get_face_neighbor_inside (tree_class, child, neigh, iface, &face_num);
    scheme->element_get_face_neighbor_inside (tree_class, neigh, element, face_num, &check);

    EXPECT_TRUE (scheme->element_is_equal (tree_class, child, element)) << "Got a false neighbor.";
    EXPECT_ELEM_EQ (scheme, tree_class, child, element);
  }
}

int
t8_test_get_middle_child (t8_eclass_t eclass, int ilevel, t8_element_t *element, t8_element_t *child, t8_scheme *scheme,
                          const t8_eclass_t tree_class)
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
    scheme->element_get_child (tree_class, element, 0, child);
    scheme->element_copy (tree_class, child, element);
    return 3;
  case T8_ECLASS_TRIANGLE:
    return 3;
  case T8_ECLASS_HEX:
    /* There are no inner children in level one refinement. The test starts with level two, because this is the first level, inner children existing.
       The third child of level one child 4 is one of eight middle children in level two. */
    scheme->element_get_child (tree_class, element, 4, child);
    scheme->element_copy (tree_class, child, element);
    return 3;
  case T8_ECLASS_TET:
    return 3;
  case T8_ECLASS_PRISM:
    /* There are no inner children in level one refinement. The test starts with level two, because this is the first level, inner children existing.
       The last child of level one child 4 is one of eight middle children in level two. */
    scheme->element_get_child (tree_class, element, 4, child);
    scheme->element_copy (tree_class, child, element);
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

  const int num_faces = scheme->element_get_num_faces (tree_class, element);

  for (int iface = 0; iface < num_faces; iface++) {

    const int num_children = scheme->element_get_num_face_children (tree_class, element, iface);
    int *child_indices = T8_ALLOC (int, num_children);
    t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
    scheme->element_new (tree_class, num_children, children);
    scheme->element_get_children_at_face (tree_class, element, iface, children, num_children, child_indices);

    for (int jchild = 0; jchild < num_children; jchild++) {

      const int child_id = child_indices[jchild];
      const int face_contact = scheme->element_face_get_child_face (tree_class, element, iface, jchild);

      scheme->element_get_child (tree_class, element, child_id, child);
      int face_num;
      int inside = scheme->element_get_face_neighbor_inside (tree_class, child, neigh, face_contact, &face_num);

      ASSERT_EQ (inside, 0) << "Element is not outside.";

      inside = scheme->element_get_tree_face (tree_class, child, face_contact);
      ASSERT_EQ (inside, iface) << "Wrong face.";
    }
    scheme->element_destroy (tree_class, num_children, children);
    T8_FREE (children);
    T8_FREE (child_indices);
  }
}

void
t8_recursive_check_diff (t8_element_t *element, t8_element_t *child, t8_element_t *neigh, t8_scheme *scheme,
                         const t8_eclass_t tree_class, int maxlvl, int level)
{

  T8_ASSERT (level <= maxlvl && maxlvl <= scheme->get_maxlevel (tree_class) - 1);
  if (level == maxlvl) {
    return;
  }

  /* Compute the neighbors neighbor along a given face and check, if the result is the
   * original element. */
  const int num_faces = scheme->element_get_num_faces (tree_class, element);

  t8_test_face_neighbor_inside (num_faces, child, element, neigh, scheme, tree_class);

  const int num_children = scheme->element_get_num_children (tree_class, child);
  for (int ichild = 0; ichild < num_children; ichild++) {
    scheme->element_get_child (tree_class, element, ichild, child);
    t8_recursive_check_diff (child, element, neigh, scheme, tree_class, maxlvl, level + 1);
    scheme->element_get_parent (tree_class, child, element);
  }
}

/* Recursively check, if all neighbors are computed correct up to a given level. */
TEST_P (face_neigh, recursive_check_diff)
{
  int level = 1;
  const int middle_child_id = t8_test_get_middle_child (tree_class, level, element, child, scheme, tree_class);
  scheme->element_get_child (tree_class, element, middle_child_id, child);

  t8_recursive_check_diff (child, element, neigh, scheme, tree_class, maxlvl, level);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_neigh, face_neigh, AllEclasses, print_eclass);
