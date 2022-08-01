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
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <sc_functions.h>
#include <t8_cmesh/t8_cmesh_examples.h>

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
}

static int
t8_test_init_linear_id_refine_everything (t8_forest_t forest,
                                          t8_forest_t forest_from,
                                          t8_locidx_t which_tree,
                                          t8_locidx_t lelement_id,
                                          t8_eclass_scheme_c *ts,
                                          const int is_family,
                                          const int num_elements,
                                          t8_element_t *elements[])
{

  return 1;
}

static void
t8_check_uniform_forest (t8_eclass_scheme_c *ts, t8_scheme_cxx_t *scheme,
                         sc_MPI_Comm comm, int maxlvl)
{
  t8_forest_t         forest, forest_adapt;
  t8_cmesh_t          cmesh;
  t8_locidx_t         first_id, last_id, j, id;
  t8_locidx_t         first_tid, last_tid, tree_id;
  t8_element_t       *element;
  int                 i;

  cmesh = t8_cmesh_new_from_class (ts->eclass, comm);
  t8_cmesh_ref (cmesh);
  forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  t8_scheme_cxx_ref (scheme);
  for (i = 0; i < maxlvl; i++) {

    /*Get the id of the first and last tree on this process */
    first_tid = t8_forest_get_first_local_tree_id (forest);
    last_tid = first_tid + t8_forest_get_num_local_trees (forest);
    /*Iterate over trees */
    for (tree_id = first_id; tree_id < last_tid; tree_id++) {
      /*Get id of the first and last element on this tree */
      first_id = t8_forest_get_first_local_element_id (forest);
      last_id = first_id + t8_forest_get_local_num_elements (forest);
      /*Iterate over elements */
      for (j = first_id; j < last_id; j++) {
        /*Get the j-th element and check the computed linear id */
        element = t8_forest_get_element (forest, j, &tree_id);
        /*Get the ID of the element at current level */
        id = ts->t8_element_get_linear_id (element, i);
        SC_CHECK_ABORT (id == j, "Wrong ID\n");
      }
    }
    t8_forest_init (&forest_adapt);
    t8_forest_set_level (forest_adapt, i + 1);
    t8_forest_set_adapt (forest_adapt, forest,
                         t8_test_init_linear_id_refine_everything, 0);
    t8_forest_commit (forest_adapt);
    forest = forest_adapt;
    t8_debugf ("Done with eclass %s at level %i\n",
               t8_eclass_to_string[ts->eclass], i);
  }
  t8_cmesh_unref (&cmesh);
  t8_forest_unref (&forest_adapt);
  t8_debugf ("Done uniform forest\n");
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
    ts->t8_element_count_leafs (element, level);
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
  t8_scheme_cxx_t    *scheme;
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
    t8_productionf ("Testing linear_id for eclass %s\n",
                    t8_eclass_to_string[eclassi]);

    eclass = (t8_eclass_t) eclassi;
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];
    t8_check_uniform_forest (ts, scheme, sc_MPI_COMM_WORLD, maxlvl);
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
  t8_init (SC_LP_DEFAULT);

  t8_check_linear_id (maxlvl);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
