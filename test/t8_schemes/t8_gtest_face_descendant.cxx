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
  int                 num_faces = ts->t8_element_num_faces (elem);

  /* Calculate the type (with eclass?, seperate function?) */
  int                 type;

  for (int ilevel = level + 1; ilevel < maxlvl; ilevel++) {
    t8_element_t       *children[];
    int                 num_children = ts->t8_element_num_children (elem);
    ts->t8_element_children (elem, num_children, children);

    for (int jface = 0; jface < num_faces; jface++) {

      /* Compute the child-id of the first-descendent */
      ts->t8_element_children_at_face (elem, jface, children, num_children,
                                       child_indices);
      int                 child_id = ts->t8_element_child_id (children[0]);

      /* Manually computing the first-descendant */
      ts->t8_element_first_descendant (elem, test, ilevel);

      ts->t8_element_first_descendant_face (elem, jface, tmp, ilevel);
      ASSERT_FALSE (ts->t8_element_compare (test,
                                            tmp)) <<
        "Wrong first descendant face\n";

      /* Computing the child-id of the last descendant */
      /* Difference between types */

      /* ts->t8_element_last_descendant (elem, test, ilevel);  

         ts->t8_element_last_descendant_face (elem, jface, tmp, ilevel);
         ASSERT_FALSE (ts->t8_element_compare (test, tmp)) << "Wrong last descendant face\n";
       */
    }
  }
}

/*Recursivly check the first and last descendant along a face*/
void
t8_recursive_face_desendant (t8_element_t *elem, t8_element_t *test,
                             t8_element_t *tmp, t8_element_t *child,
                             t8_eclass_scheme_c *ts, int maxlvl)
{

}

TEST_P (class_descendant, t8_check_face_desc)
{
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
INSTANTIATE_TEST_SUITE_P (t8_gtest_element_face_descendant, class_descendant,testing::Range(T8_ECLASS_LINE, T8_ECLASS_COUNT));
/* *INDENT-ON* */
