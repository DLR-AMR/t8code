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
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.h>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.h>

/* *INDENT-OFF* */
class class_pyra_descendant:public testing::TestWithParam <t8_eclass_t > {
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
t8_linear_face_descendant (t8_element_t *elem, t8_element_t *tmp,
                           t8_element_t *test, t8_eclass_scheme_c *ts,
                           int maxlvl)
{
  int                 level = ts->t8_element_level (elem);
  int                 type = ((t8_dpyramid_t *) elem)->pyramid.type;
  int                 child_id;
  int                 num_faces;
  
  if (type < 6) {
    num_faces = T8_DTET_FACES;
  }
  else {
    num_faces = T8_DPYRAMID_FACES;
  }
  for (int ilevel = level + 1; ilevel < maxlvl; ilevel++) {
    for (int jface = 0; jface < num_faces; jface++) {
      /* Compute the child-id of the first-descendant */
      if (type >= 6) {
        child_id = t8_dpyramid_type_face_to_children_at_face[type - 6][jface][0];
      }
      else {
        child_id = t8_dtet_face_corner[jface][0];
        child_id = t8_dtet_parenttype_beyid_to_Iloc[type][child_id];
      }
      /* Manually computing the first_descendant */
      ts->t8_element_copy (elem, tmp);
      for (int klevel = level; klevel < ilevel; klevel++) {
        ts->t8_element_child (tmp, child_id, test);
        ts->t8_element_copy (test, tmp);
      }

      ts->t8_element_first_descendant_face (elem, jface, tmp, ilevel);
      EXPECT_TRUE(!ts->t8_element_compare (test, tmp));

      /* Computing the child-id of the last descendant */
      if (type >= 6) {
        child_id = t8_dpyramid_type_face_to_children_at_face[type - 6][jface][3];
      }
      else {
        child_id = SC_MAX (t8_dtet_face_corner[jface][1], t8_dtet_face_corner[jface][2]);
        child_id = t8_dtet_parenttype_beyid_to_Iloc[type][child_id];
      }
      /* Manually computing the last_descendant */
      ts->t8_element_copy (elem, tmp);
      for (int klevel = level; klevel < ilevel; klevel++) {
        ts->t8_element_child (tmp, child_id, test);
        ts->t8_element_copy (test, tmp);
      }
      ts->t8_element_last_descendant_face (elem, jface, tmp, ilevel);

      EXPECT_TRUE(!ts->t8_element_compare (test, tmp));
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
  int                 child_id;
  int                 type;
  int                 num_faces;
  int                 num_children;
  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  if (level == maxlvl)
    return;
  /* Get the number of faces of the current element */
  if (ts->t8_element_shape (elem) == T8_ECLASS_PYRAMID) {
    num_faces = T8_DPYRAMID_FACES;
  }
  else {
    num_faces = T8_DTET_FACES;
  }
  type = ((t8_dpyramid_t *) elem)->pyramid.type;
  /* Check face descendants from current-level + 1 to given maximum level */
  for (int ilevel = level + 1; ilevel < maxlvl; ilevel++) {
    for (int jface = 0; jface < num_faces; jface++) {
      /* Get the id of the child at the face */
      if (type >= 6) {
        child_id = t8_dpyramid_type_face_to_children_at_face[type - 6][jface][0];
      }
      else {
        child_id = t8_dtet_face_corner[jface][0];
        child_id = t8_dtet_parenttype_beyid_to_Iloc[type][child_id];
      }

      ts->t8_element_copy (elem, tmp);
      /* The child-id of the descendant is always the same over the level. */
      for (int klevel = level; klevel < ilevel; klevel++) {
        ts->t8_element_child (tmp, child_id, test);
        ts->t8_element_copy (test, tmp);
      }

      ts->t8_element_first_descendant_face (elem, jface, tmp, ilevel);
      EXPECT_TRUE(!ts->t8_element_compare (test, tmp));

      /* Analogously check the last facedescendant */
      if (type >= 6) {
        child_id = t8_dpyramid_type_face_to_children_at_face[type - 6][jface][3];
      }
      else {
        child_id =
          SC_MAX (t8_dtet_face_corner[jface][1], t8_dtet_face_corner[jface][2]);
        child_id = t8_dtet_parenttype_beyid_to_Iloc[type][child_id];
      }
      ts->t8_element_copy (elem, tmp);
      for (int klevel = level; klevel < ilevel; klevel++) {
        ts->t8_element_child (tmp, child_id, test);
        ts->t8_element_copy (test, tmp);
      }
      ts->t8_element_last_descendant_face (elem, jface, tmp, ilevel);
      EXPECT_TRUE(!ts->t8_element_compare (test, tmp));
    }
  }

  num_children = ts->t8_element_num_children (elem);
  for (int ichild = 0; ichild < num_children; ichild++) {
    ts->t8_element_child (elem, ichild, child);
    t8_recursive_face_desendant (child, test, tmp, elem, ts, maxlvl);
    ts->t8_element_parent (child, elem);
  }
}

TEST_P(class_pyra_descendant,  t8_check_face_desc ){
#ifdef T8_ENABLE_DEBUG
    const int           maxlvl = 3;
#else
    const int           maxlvl = 4;
#endif
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
INSTANTIATE_TEST_SUITE_P (t8_gtest_element_pyra_face_descendant, class_pyra_descendant,testing::Range(T8_ECLASS_PYRAMID, T8_ECLASS_COUNT));
/* *INDENT-ON* */