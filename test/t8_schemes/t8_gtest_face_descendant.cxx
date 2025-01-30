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
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_gtest_dfs_base.hxx"

/* compute the first/last descendant by iteratively taking the first/last child at each level*/
static void
t8_test_manual_first_last_face_descendant (const t8_scheme *scheme, const t8_element_t *element,
                                           const t8_eclass_t tree_class, const int iface, const int desc_level,
                                           const int last, t8_element_t *face_desc)
{
  const int num_children_at_face = scheme->element_get_num_face_children (tree_class, element, iface);

  int *child_indices = T8_ALLOC (int, num_children_at_face);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children_at_face);
  scheme->element_new (tree_class, num_children_at_face, children);

  scheme->element_copy (tree_class, element, face_desc);
  const int level = scheme->element_get_level (tree_class, element);
  for (int ilevel = level; ilevel < desc_level; ilevel++) {
    EXPECT_EQ (scheme->element_get_num_face_children (tree_class, element, iface), num_children_at_face);
    /* Compute child_id of the test_child_id-th child. */
    scheme->element_get_children_at_face (tree_class, face_desc, iface, children, num_children_at_face, child_indices);

    /* chose correct face_id dependent on if we want first or last face desc.*/
    const int face_child_id = last ? num_children_at_face - 1 : 0;

    const int child_id = child_indices[face_child_id];

    scheme->element_get_child (tree_class, face_desc, child_id, face_desc);
  }
  scheme->element_destroy (tree_class, num_children_at_face, children);
  T8_FREE (children);
  T8_FREE (child_indices);
}

class class_descendant: public TestDFS {
  void
  check_element () override
  {
    /* Check the linear first and last descendants of an element along all faces. 
     * For the test the descendants are computed manually by t8_test_manual_first_last_face_descendant and 
     * by the scheme implementation t8_element_first_descendant for the first descendants over the levels.
     */

    const int level = scheme->element_get_level (eclass, element);
    const int num_faces = scheme->element_get_num_faces (eclass, element);

    /* Testing the linear first descendant. */
    for (int ilevel = level + 1; ilevel < max_test_lvl; ilevel++) {
      for (int jface = 0; jface < num_faces; jface++) {

        t8_test_manual_first_last_face_descendant (scheme, element, eclass, jface, ilevel, 0, manual_face_desc);
        scheme->element_get_first_descendant_face (eclass, element, jface, scheme_face_desc, ilevel);
        /* Compare the manually computed child with the result of t8_element_first_descendant_face. */
        EXPECT_ELEM_EQ (scheme, eclass, scheme_face_desc, manual_face_desc);

        t8_test_manual_first_last_face_descendant (scheme, element, eclass, jface, ilevel, 1, manual_face_desc);
        scheme->element_get_last_descendant_face (eclass, element, jface, scheme_face_desc, ilevel);
        /* Compare the manually computed child with the result of t8_element_last_descendant_face. */
        EXPECT_ELEM_EQ (scheme, eclass, scheme_face_desc, manual_face_desc);
      }
    }
  }

 protected:
  void
  SetUp () override
  {
    dfs_test_setup ();
    max_test_lvl = scheme->get_maxlevel (eclass);
    scheme->element_new (eclass, 1, &manual_face_desc);
    scheme->element_new (eclass, 1, &scheme_face_desc);
  }
  void
  TearDown () override
  {
    scheme->element_destroy (eclass, 1, &manual_face_desc);
    scheme->element_destroy (eclass, 1, &scheme_face_desc);
    dfs_test_teardown ();
  }
  int max_test_lvl;
  t8_element_t *manual_face_desc;
  t8_element_t *scheme_face_desc;
};

TEST_P (class_descendant, t8_check_face_desc)
{

#ifdef T8_ENABLE_LESS_TESTS
  const int maxlvl = 3;
#else
  const int maxlvl = 5;
#endif

  check_recursive_dfs_to_max_lvl (maxlvl);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_face_descendant, class_descendant, AllSchemes);
