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

#include <t8_eclass.h>
#include <t8_schemes/t8_default_cxx.hxx>


/* Build every child for the parent element and check wether the computed parent is
 * again the input element. Then call this function for every child, until maxlevel is
 * reached.
 */
static void
t8_recursive_child_find_parent (t8_element_t * element,
                                t8_eclass_scheme_c * ts, int level,
                                int maxlvl)
{
  T8_ASSERT (level <= maxlvl
             || maxlvl < ts->t8_element_num_children (element));
  int                 num_children, i;
  t8_element_t       *child, *test_parent;
  /*Get number of children */
  num_children = ts->t8_element_num_children (element);
  /*Get child and test_parent, to check if test_parent = parent of child */
  ts->t8_element_new (1, &child);
  ts->t8_element_new (1, &test_parent);
  if (level == maxlvl)
    return;
  for (i = 0; i < num_children; i++) {
    /*Compute child i */
    ts->t8_element_child (element, i, child);
    /*Compute parent of child */
    ts->t8_element_parent (child, test_parent);
    /*If its equal, call child_find_parent, to check if parent-child relation
     * is correct in next level until maxlvl is reached*/
    if (ts->t8_element_compare (element, test_parent)) {
      SC_ABORT ("Computed child_parent is not the parent");
    }
    else {
      t8_recursive_child_find_parent (child, ts, level + 1, maxlvl);
    }
  }
}

static void
t8_compute_child_find_parent (int maxlvl)
{
  t8_element_t       *element;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  t8_eclass_t         eclass;
  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_PYRAMID; eclassi++) {
    eclass = (t8_eclass_t) eclassi;
    /*Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
    /*Get element and initialize it */
    ts->t8_element_new (1, &element);
    ts->t8_element_set_linear_id (element, 0, 0);
    /*Check for correct parent-child relation */
    t8_recursive_child_find_parent (element, ts, 0, maxlvl);
    printf ("%s: Success\n", t8_eclass_to_string[eclass]);
    /*Destroy element */
    ts->t8_element_destroy (1, &element);

  }
  /*Destroy scheme */
  t8_scheme_cxx_unref (&scheme);
}

int
main ()
{
  t8_compute_child_find_parent (4);
  return 0;
}
