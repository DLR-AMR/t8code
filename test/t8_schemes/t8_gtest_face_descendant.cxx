/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/*TODO: Delete GTEST_SKIP in TEST_P t8_check_face_desc if t8_prism_children_at_face is fixed. */

/* *INDENT-OFF* */
class class_descendant:public testing::TestWithParam <t8_eclass_t > {
protected:
  void SetUp () override {
    eclass = GetParam();

    scheme = t8_scheme_new_default_cxx ();
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
  }
  void TearDown () override {
    /* Destroy scheme */
    t8_scheme_cxx_unref (&scheme);
  }
  t8_eclass_t         eclass;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;

};
/* *INDENT-ON* */

void
t8_face_descendant_test_child (t8_element_t *tmp, t8_element_t *test,
                               t8_element_t *elem, t8_eclass_scheme_c *ts,
                               int face, int level, int ilevel,
                               int num_children, int test_child_id)
{

  int                *child_indices = T8_ALLOC (int, num_children);

  t8_element_t      **children = T8_ALLOC (t8_element_t *, num_children);
  ts->t8_element_new (num_children, children);

  /* Computing the child-id of the last descendant */
  ts->t8_element_children_at_face (tmp, face, children, num_children,
                                   child_indices);
  int                 child_id = child_indices[test_child_id];
  ts->t8_element_copy (elem, tmp);
  for (int klevel = level; klevel < ilevel; klevel++) {
    ts->t8_element_child (tmp, child_id, tmp);
  }
  ts->t8_element_destroy (num_children, children);
  T8_FREE (children);
  T8_FREE (child_indices);
}

void
t8_linear_face_descendant (t8_element_t *elem, t8_element_t *tmp,
                           t8_element_t *test, t8_eclass_scheme_c *ts,
                           int maxlvl)
{

  int                 level = ts->t8_element_level (elem);
  int                 num_faces;

  num_faces = ts->t8_element_num_faces (elem);

  ts->t8_element_copy (elem, tmp);
  for (int ilevel = level + 1; ilevel < maxlvl; ilevel++) {
    for (int jface = 0; jface < num_faces; jface++) {

      int                 num_children =
        ts->t8_element_num_face_children (tmp, jface);
      t8_face_descendant_test_child (tmp, test, elem, ts, jface, level,
                                     ilevel, num_children, 0);

      ts->t8_element_first_descendant_face (elem, jface, tmp, ilevel);
      ASSERT_FALSE (ts->t8_element_compare (test,
                                            tmp)) <<
        "Wrong first descendant face\n";
    }
  }

  ts->t8_element_copy (elem, tmp);
  for (int ilevel = level + 1; ilevel < maxlvl; ilevel++) {
    for (int jface = 0; jface < num_faces; jface++) {

      int                 num_children =
        ts->t8_element_num_face_children (tmp, jface);
      t8_face_descendant_test_child (tmp, test, elem, ts, jface, level,
                                     ilevel, num_children, num_children - 1);

      ts->t8_element_last_descendant_face (elem, jface, tmp, ilevel);
      ASSERT_FALSE (ts->t8_element_compare (test,
                                            tmp)) <<
        "Wrong last descendant face\n";
    }
  }
}

/*Recursivly check the first and last descendant along a face*/
void
t8_recursive_face_desendant (t8_element_t *elem, t8_element_t *test,
                             t8_element_t *tmp, t8_element_t *child,
                             t8_eclass_scheme_c *ts, int maxlvl)
{
  int                 level = ts->t8_element_level (elem);
  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  if (level == maxlvl)
    return;

  t8_linear_face_descendant (elem, tmp, test, ts, maxlvl);

  int                 num_children = ts->t8_element_num_children (elem);
  for (int ichild = 0; ichild < num_children; ichild++) {
    ts->t8_element_child (elem, ichild, child);
    t8_recursive_face_desendant (child, test, tmp, elem, ts, maxlvl);
    ts->t8_element_parent (child, elem);
  }
}

TEST_P (class_descendant, t8_check_face_desc)
{

#ifdef T8_ENABLE_DEBUG
  const int           maxlvl = 3;
#else
  const int           maxlvl = 4;
#endif
  if ((int) eclass == (int) T8_ECLASS_PRISM) {
    GTEST_SKIP ();
  }
  t8_element_t       *element;
  t8_element_t       *child;
  t8_element_t       *test;
  t8_element_t       *tmp;

  /* Get element and initialize it */
  ts->t8_element_new (1, &element);
  ts->t8_element_new (1, &child);
  ts->t8_element_new (1, &test);
  ts->t8_element_new (1, &tmp);

  ts->t8_element_set_linear_id (element, 0, 0);

  /* Check for correct parent-child relation */
  t8_linear_face_descendant (element, child, test, ts, maxlvl);
  t8_recursive_face_desendant (element, test, tmp, child, ts, maxlvl);

  /* Destroy element */
  ts->t8_element_destroy (1, &element);
  ts->t8_element_destroy (1, &child);
  ts->t8_element_destroy (1, &test);
  ts->t8_element_destroy (1, &tmp);
}

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_element_face_descendant, class_descendant,testing::Range(T8_ECLASS_VERTEX, T8_ECLASS_COUNT));
/* *INDENT-ON* */
