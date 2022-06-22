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
#include <sc_functions.h>

/*recursively compute all elements and check their id*/
static void
t8_recursive_linear_id (t8_element_t *element, t8_element_t *child,
                        t8_element_t *test, t8_eclass_scheme_c *ts,
                        int maxlvl, uint64_t *id)
{
  int                 level = ts->t8_element_level (element);
  int                 num_children, i;
  num_children = ts->t8_element_num_children (element);
  if (level == maxlvl - 1) {
    for (i = 0; i < num_children; i++) {
      ts->t8_element_child (element, i, child);
      ts->t8_element_set_linear_id (test, maxlvl, *id);

      SC_CHECK_ABORT (!ts->t8_element_compare (child, test),
                      "Wrong element\n");
      SC_CHECK_ABORT (ts->t8_element_get_linear_id (test, maxlvl) == (*id),
                      "Wrong linear id\n");
      (*id)++;
    }
  }
  else {
    for (i = 0; i < num_children; i++) {
      ts->t8_element_child (element, i, child);
      t8_recursive_linear_id (child, element, test, ts, maxlvl, id);
      ts->t8_element_parent (child, element);
    }
  }
}

int
t8_num_descendants (t8_element_t *element, int level, t8_eclass_scheme_c *ts)
{
  int                 level_diff = level - ts->t8_element_level (element);
  int                 shape = ts->t8_element_shape (element);
  switch (shape) {
  case T8_ECLASS_VERTEX:
    return 1;
  case T8_ECLASS_LINE:
    return sc_intpow (2, level_diff);
  case T8_ECLASS_QUAD:
    return sc_intpow (4, level_diff);
  case T8_ECLASS_TRIANGLE:
    return sc_intpow (4, level_diff);
  case T8_ECLASS_HEX:
    return sc_intpow (8, level_diff);
  case T8_ECLASS_TET:
    return sc_intpow (8, level_diff);
  case T8_ECLASS_PRISM:
    return sc_intpow (8, level_diff);
  case T8_ECLASS_PYRAMID:
    return 2 * sc_intpow (8, level_diff) - sc_intpow (6, level_diff);
  default:
    SC_ABORT ("Class Non-existent");
  }
}

/*Check, if all descendants of an element at level maxlvl have the same id on
 * the level of the input element as the input element*/
static void
t8_id_at_other_lvl_check (t8_element_t *element, t8_element * child,
                          t8_eclass_scheme_c *ts, int maxlvl)
{
  int                 level = ts->t8_element_level (element);
  t8_linearidx_t      current_id =
    ts->t8_element_get_linear_id (element, level);
  t8_linearidx_t      num_descendants =
    t8_num_descendants (element, maxlvl, ts);
  t8_linearidx_t      id_at_lvl =
    ts->t8_element_get_linear_id (element, maxlvl);
  t8_linearidx_t      i;
  t8_linearidx_t      id;

  for (i = 0; i < num_descendants; i++) {
    ts->t8_element_set_linear_id (child, maxlvl, id_at_lvl + i);
    id = ts->t8_element_get_linear_id (child, level);
    SC_CHECK_ABORT (id == current_id, "Wrong id");
  }

}

static void
t8_check_linear_id (const int maxlvl)
{
  t8_element_t       *element, *child, *test;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  int                 level;
  int                 num_desc;
  int                 i;
  int                 j;
  t8_linearidx_t      id;
  t8_eclass_t         eclass;

  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_COUNT; eclassi++) {
    eclass = (t8_eclass_t) eclassi;
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
    /* Get element and initialize it */
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test);

    ts->t8_element_set_linear_id (element, 0, 0);
    /* Check for correct parent-child relation */
    for (level = 1; level <= maxlvl; level++) {
      id = 0;
      t8_recursive_linear_id (element, child, test, ts, level, &id);
    }
    if (eclassi == T8_ECLASS_PYRAMID) {
      ts->t8_element_set_linear_id (element, 0, 0);
      for (j = 0; j < maxlvl - 4; j++) {
        num_desc = t8_num_descendants (element, j, ts);
        for (i = 0; i < num_desc; i++) {
          ts->t8_element_set_linear_id (child, j, i);
          t8_id_at_other_lvl_check (child, test, ts, maxlvl);
        }
      }
    }
    /* Destroy element */
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &test);

  }
  /* Destroy scheme */
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

  t8_check_linear_id (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
