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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/*This program tests the descendant function of an element.*/

/* Test recursively if the first and last descendant of an element is
 * computed correctly. Only the descendant of elem->level + 1 is tested.*/
static void
t8_recursive_descendant (t8_element_t *elem, t8_element_t *desc,
                         t8_element_t *test, t8_eclass_scheme_c *ts,
                         int maxlvl)
{
  int                 num_children = ts->t8_element_num_children (elem);
  int                 level = ts->t8_element_level (elem);
  int                 i;
  for (i = 0; i < num_children; i++) {
    ts->t8_element_child (elem, i, desc);
    /*first child == first descendant */
    if (i == 0) {
      ts->t8_element_first_descendant (elem, test, level + 1);
      SC_CHECK_ABORT (!ts->t8_element_compare (desc, test),
                      "Wrong first descendant\n");
    }
    /*last child == last descendant */
    else if (i == num_children - 1) {
      ts->t8_element_last_descendant (elem, test, level + 1);
      SC_CHECK_ABORT (!ts->t8_element_compare (desc, test),
                      "Wrong last descendant\n");
    }
    else if (level > maxlvl) {
      t8_recursive_descendant (desc, elem, test, ts, maxlvl);
      ts->t8_element_parent (desc, elem);
    }
  }
}

/* Test, if the first descendant of an element is computed correctly over a range
 * of levels*/
static void
t8_deep_first_descendant (t8_element_t *elem, t8_element_t *desc,
                          t8_element_t *test, t8_eclass_scheme_c *ts,
                          int level)
{
  int                 i, elem_level = ts->t8_element_level (elem);
  ts->t8_element_copy (elem, test);
  for (i = elem_level; i < level; i++) {
    ts->t8_element_child (test, 0, desc);
    ts->t8_element_copy (desc, test);
  }
  ts->t8_element_first_descendant (elem, test, level);
  SC_CHECK_ABORT (!ts->t8_element_compare (desc, test),
                  "Wrong deep first descendant\n");
}

/* Test, if the last descendant of an element is computed correctly over a range
 * of levels*/
static void
t8_deep_last_descendant (t8_element_t *elem, t8_element_t *desc,
                         t8_element_t *test, t8_eclass_scheme_c *ts,
                         int level)
{
  int                 num_children, i;
  ts->t8_element_copy (elem, test);
  /*Compute the correct element */
  for (i = ts->t8_element_level (elem); i < level; i++) {
    num_children = ts->t8_element_num_children (test);
    ts->t8_element_child (test, num_children - 1, desc);
    ts->t8_element_copy (desc, test);
  }
  /*Check for equality */
  ts->t8_element_last_descendant (elem, test, level);
  SC_CHECK_ABORT (!ts->t8_element_compare (desc, test),
                  "Wrong deep last descendant\n");
}

/* Test if the first and last descendant of an element are computed correctly.
 * The level between the element and the descendant is larger or equal to one */
static void
t8_large_step_descendant (t8_element_t *elem, t8_element_t *desc,
                          t8_element_t *test, t8_eclass_scheme_c *ts,
                          int maxlvl)
{
  int                 num_children = ts->t8_element_num_children (elem);
  int                 i, j;
  for (i = ts->t8_element_level (elem); i < maxlvl; i++) {
    /*Use these functions to perform the actual test */
    t8_deep_first_descendant (elem, desc, test, ts, maxlvl);
    t8_deep_last_descendant (elem, desc, test, ts, maxlvl);
    for (j = 0; j < num_children; j++) {
      ts->t8_element_child (elem, j, desc);
      t8_large_step_descendant (desc, elem, test, ts, maxlvl);
      ts->t8_element_parent (desc, elem);
    }
  }
}

static void
t8_check_descendant (const int maxlvl)
{
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  t8_eclass_t         eclass;
  t8_element_t       *elem, *desc, *test;
  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_COUNT; eclassi++) {
    t8_productionf ("Checking descendants for eclass %s\n",
                    t8_eclass_to_string[eclassi]);
    eclass = (t8_eclass_t) eclassi;
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &elem);
    ts->t8_element_new (1, &desc);
    ts->t8_element_new (1, &test);
    ts->t8_element_set_linear_id (elem, 0, 0);

    t8_recursive_descendant (elem, desc, test, ts, maxlvl);
    t8_deep_first_descendant (elem, desc, test, ts,
                              ts->t8_element_maxlevel ());
    t8_deep_last_descendant (elem, desc, test, ts,
                             ts->t8_element_maxlevel ());
    t8_large_step_descendant (elem, desc, test, ts, maxlvl);

    ts->t8_element_destroy (1, &elem);
    ts->t8_element_destroy (1, &desc);
    ts->t8_element_destroy (1, &test);
  }
  t8_scheme_cxx_unref (&scheme);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
#ifdef T8_ENABLE_DEBUG
  const int           maxlvl = 3;
#else
  const int           maxlvl = 4;
#endif

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_check_descendant (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
