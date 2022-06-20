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

#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

/* Build every child for the parent element and check wether the computed parent is
 * again the input element. Then call this function for every child, until maxlevel is
 * reached.
 */
static void
t8_recursive_child_find_parent (t8_element_t *element, t8_element_t *child,
                                t8_element_t *test_parent,
                                t8_eclass_scheme_c *ts, int level,
                                const int maxlvl)
{
  T8_ASSERT (level <= maxlvl && maxlvl <= ts->t8_element_maxlevel () - 1);
  int                 num_children, i;
  /* Get number of children */
  num_children = ts->t8_element_num_children (element);
  /* Get child and test_parent, to check if test_parent = parent of child */
  if (level == maxlvl)
    return;
  for (i = 0; i < num_children; i++) {
    /* Compute child i */
    ts->t8_element_child (element, i, child);
    /* Compute parent of child */
    ts->t8_element_parent (child, test_parent);
    /* If its equal, call child_find_parent, to check if parent-child relation
     * is correct in next level until maxlvl is reached*/
    SC_CHECK_ABORT (!ts->t8_element_compare (element, test_parent),
                    "Computed child_parent is not the parent");

    t8_recursive_child_find_parent (child, element, test_parent, ts,
                                    level + 1, maxlvl);
    /* After the check we know the parent-function is correct for this child.
     * Therefore we can use it to recompute the element*/
    ts->t8_element_parent (child, element);
  }
}

static void
t8_compute_child_find_parent (const int maxlvl)
{
  t8_element_t       *element, *child, *test_parent;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  t8_eclass_t         eclass;
  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_COUNT; eclassi++) {
    t8_productionf ("Testing child find parent for eclass %s\n",
                    t8_eclass_to_string[eclassi]);
    eclass = (t8_eclass_t) eclassi;
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
    /* Get element and initialize it */
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test_parent);

    ts->t8_element_set_linear_id (element, 0, 0);
    /* Check for correct parent-child relation */
    t8_recursive_child_find_parent (element, child, test_parent, ts, 0,
                                    maxlvl);
    /* Destroy element */
    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &test_parent);

  }
  /* Destroy scheme */
  t8_scheme_cxx_unref (&scheme);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
#ifdef T8_ENABLE_LESS_TESTS
  const int           maxlvl = 4;
#else
  const int           maxlvl = 9;
#endif

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_compute_child_find_parent (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
