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

/* Check the computation of the successor recursively. Iterate through the elements
 * via DFS. On the given maxlvl-1 the children are computeted iteratively. For
 * each child, the successor is checked.
 */
static void
t8_recursive_successor (t8_element_t *element, t8_element_t *successor,
                        t8_element_t *child,
                        t8_element_t *last, t8_eclass_scheme_c *ts,
                        const int maxlvl)
{
  int                 level = ts->t8_element_level (element);
  T8_ASSERT (ts->t8_element_level (element) <= maxlvl
             && maxlvl <= ts->t8_element_maxlevel () - 1);
  int                 num_children, i;
  num_children = ts->t8_element_num_children (element);
  if (level == maxlvl - 1) {
    /* Check, if the successor of the last recursion is the first child of
     * of this element.
     */
    ts->t8_element_child (element, 0, child);
    SC_CHECK_ABORT (!ts->t8_element_compare (child, successor),
                    "Wrong Successor, Case1\n");
    num_children = ts->t8_element_num_children (element);
    /*Check if the successor in this element is computed correctly */
    for (i = 1; i < num_children; i++) {
      ts->t8_element_successor (child, successor, maxlvl);
      ts->t8_element_child (element, i, child);
      SC_CHECK_ABORT (!ts->t8_element_compare (child, successor),
                      "Wrong Succesor, Case2\n");
    }
    /*If the iterator is the last element, the test can finish */
    if (!ts->t8_element_compare (last, child)) {
      return;
    }
    /*Compute the next successor / "jump" out of the current element */
    else {
      ts->t8_element_successor (child, successor, maxlvl);
    }
  }
  else {
    /*DFS run through the elements */
    num_children = ts->t8_element_num_children (element);
    for (i = 0; i < num_children; i++) {
      ts->t8_element_child (element, i, child);
      t8_recursive_successor (child, successor, element, last, ts, maxlvl);
      ts->t8_element_parent (child, element);
    }
  }
}

/* Check the computation of the successor at the maximum level of the element.
 * Given the element of level maxlevel-2  with linear id 0 all children of
 * maximum level are computed. The successor runs through all these children.
 */
static void
t8_deep_successor (t8_element_t *element, t8_element_t *successor,
                   t8_element_t *child, t8_eclass_scheme_c *ts)
{
  int                 maxlvl = ts->t8_element_maxlevel ();
  int                 i, j, num_children, num_children_child;
  num_children = ts->t8_element_num_children (element);
  for (i = 0; i < num_children; i++) {
    ts->t8_element_child (element, i, child);
    /*Go to the children at maximum level */
    num_children_child = ts->t8_element_num_children (child);
    for (j = 0; j < num_children_child; j++) {
      ts->t8_element_child (child, j, element);
      /*Check the computation of the successor */
      SC_CHECK_ABORT (!ts->t8_element_compare (element, successor)
                      , "Wrong Successor at Maxlvl\n");
      /*Compute the next successor */
      ts->t8_element_successor (successor, successor, maxlvl);
    }
    ts->t8_element_parent (child, element);
  }
}

static void
t8_compute_successor (const int level)
{
  t8_element_t       *element, *successor, *child, *last;
  t8_scheme_cxx      *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi, i;
  t8_eclass_t         eclass;
  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_LINE; eclassi < T8_ECLASS_COUNT; eclassi++) {
    eclass = (t8_eclass_t) eclassi;
    ts = scheme->eclass_schemes[eclass];
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &successor);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &last);

    ts->t8_element_set_linear_id (element, 0, 0);

    /*Test at lower level */
    for (i = 1; i <= level; i++) {
      ts->t8_element_set_linear_id (successor, i, 0);
      ts->t8_element_last_descendant (element, last, i);
      t8_recursive_successor (element, successor, child, last, ts, i);
    }
    /*Test at Maxlevel */
    ts->t8_element_set_linear_id (element, ts->t8_element_maxlevel () - 2, 0);
    ts->t8_element_set_linear_id (successor, ts->t8_element_maxlevel (), 0);
    t8_deep_successor (element, successor, last, ts);

    ts->t8_element_destroy (1, &element);
    ts->t8_element_destroy (1, &successor);
    ts->t8_element_destroy (1, &child);
    ts->t8_element_destroy (1, &last);
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

  t8_compute_successor (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
