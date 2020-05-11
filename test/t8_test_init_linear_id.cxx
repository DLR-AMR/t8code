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
#include <t8_forest.h>
#include <t8_forest/t8_forest_private.h>


static void
t8_serial_linear_id(t8_element_t * element, t8_eclass_scheme_c * ts,
                    int maxlvl)
{
    t8_gloidx_t     num_leafs;
    t8_linearidx_t  id, check_id;
    int i;
    /*Check every level up to maxlvl*/
    for(i = 0; i <= maxlvl; i++)
    {
        /*Get the number of elements with level i*/
        num_leafs = ts->t8_element_count_leafs_from_root(i);
        T8_ASSERT((t8_linearidx_t)(num_leafs) == num_leafs);
        for(id = 0; id < num_leafs; id++)
        {
            /*Construct element #id and check, if the computed id is the same*/
            ts->t8_element_set_linear_id(element, i, id);
            check_id = ts->t8_element_get_linear_id(element, i);
            SC_CHECK_ABORT(check_id == id, "Wrong ID\n");
        }
    }
}

static void
t8_check_uniform_forest(t8_eclass_scheme_c * ts, t8_scheme_cxx_t * scheme,
                        sc_MPI_Comm comm, int maxlvl)
{
    t8_forest_t         forest;
    t8_cmesh_t          cmesh;
    t8_locidx_t         first_id, last_id, j, id;
    t8_locidx_t         first_tid, last_tid, tree_id;
    t8_element_t        *element;
    int                 i;

    cmesh = t8_cmesh_new_from_class(ts->eclass, comm);
    t8_cmesh_ref(cmesh);
    for(i = 0; i < maxlvl; i++)
    {
        forest = t8_forest_new_uniform(cmesh, scheme, i, 0, comm);
        /*Reuse the forest*/
        t8_forest_ref(forest);
        /*Get the id of the first and last tree on this process*/
        first_tid = t8_forest_get_first_local_tree_id(forest);
        last_tid = first_tid + t8_forest_get_num_local_trees(forest);
        /*Iterate over trees*/
        for(tree_id = first_id; tree_id < last_tid; tree_id++)
        {
            /*Get id of the first and last element on this tree*/
            first_id = t8_forest_get_first_local_element_id(forest);
            last_id = first_id + t8_forest_get_num_element(forest);
            /*Iterate over elements*/
            for(j = first_id; j < last_id; j++){
                /*Get the j-th element and check the computed linear id*/
                element = t8_forest_get_element(forest, j, &tree_id);
                id = ts->t8_element_get_linear_id(element, i);
                SC_CHECK_ABORT(id == j, "Wrong ID\n");
            }
        }
        t8_forest_unref(&forest);
        t8_debugf("Done with eclass %s at level %i\n",
                  t8_eclass_to_string[ts->eclass], i);
    }
    t8_cmesh_unref(&cmesh);
    t8_cmesh_destroy(&cmesh);
}

static void
t8_check_linear_id (const int maxlvl)
{
  t8_element_t       *element, *child, *test;
  t8_scheme_cxx_t    *scheme;
  t8_eclass_scheme_c *ts;
  int                 eclassi;
  t8_eclass_t         eclass;

  scheme = t8_scheme_new_default_cxx ();
  for (eclassi = T8_ECLASS_ZERO; eclassi < T8_ECLASS_PYRAMID; eclassi++) {
    /* TODO: Include pyramids as soon as they are supported. */
    eclass = (t8_eclass_t) eclassi;
    /* Get scheme for eclass */
    ts = scheme->eclass_schemes[eclass];

    /* Get element and initialize it */
    ts->t8_element_new (1, &element);
    ts->t8_element_new (1, &child);
    ts->t8_element_new (1, &test);

    ts->t8_element_set_linear_id (element, 0, 0);
    /* Check for correct parent-child relation */
    t8_serial_linear_id(element, ts, maxlvl);
    t8_check_uniform_forest(ts, scheme, sc_MPI_COMM_WORLD, maxlvl);

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
  const int           maxlvl = 8;
#else
  const int           maxlvl = 9;
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
