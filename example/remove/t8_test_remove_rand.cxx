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

//#include<cstdlib>
#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "t8_cmesh/t8_cmesh_testcases.h"


static int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  if (rand()%4) {
      return -2;
  }
  return 0;
}

static int
t8_adapt_callback_coarse (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  if (is_family && rand()%10) {
    return -1;
  }
  return 0;
}

static int
t8_adapt_callback_refine (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  if (rand()%4) {
      return 1;
  }
  return 0;
}


static t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 int recursive)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  t8_forest_set_partition (forest_new, NULL, 0);
  t8_forest_commit (forest_new);

  return forest_new;
}

void
t8_test_emelemts_remove (int cmesh_id)
{
  int                 level, min_level, max_level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();

  /* Construct a cmesh */
  //cmesh = t8_test_create_cmesh (cmesh_id);
  cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 1, 0);
  //cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Compute the first level, such that no process is empty */
  min_level = t8_forest_min_nonempty_level (cmesh, scheme);

  min_level = SC_MAX (min_level, 3);
  max_level = min_level + 1;
  
  for (level = min_level; level < max_level; level++) {
    t8_cmesh_ref (cmesh);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

    for (int i = 0; i < 3; i++) {
        forest = t8_adapt_forest (forest, t8_adapt_callback_refine, 0);
        forest = t8_adapt_forest (forest, t8_adapt_callback_remove, 0);
    }
    //t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_test_remove_rand_base");

    /* coarsing until level 0 */
    // will get replaced by recursive coarseening
    for (int i = 0; i < 1*level; i++)
    {
      forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 0);
      //forest = t8_forest_new_adapt (forest, t8_adapt_callback_coarse_all, 0, 0, &adapt_data);
    }
    
    t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_test_remove_rand");
    //SC_CHECK_ABORT (t8_forest_no_overlap(forest),
    //            "The forest has overlapping elements");

    //SC_CHECK_ABORT (t8_forest_no_overlap_process_boundary (forest),
    //            "The forest has overlapping elements on process boundary");

    t8_scheme_cxx_ref (scheme);
    t8_forest_unref (&forest);
  }
  t8_scheme_cxx_unref (&scheme);
  t8_cmesh_destroy (&cmesh);
}

void
test_cmesh_emelemts_remove_all ()
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    /* This if statement is necessary to make the test work by avoiding specific cmeshes which do not work yet for this test.
     * When the issues are gone, remove the if statement. */
    if (cmesh_id != 6 && cmesh_id != 89 && (cmesh_id < 237 || cmesh_id > 256)) {
      /* Skip all t8_test_create_new_bigmesh_cmesh 
       * When issue #213 is fixed, remove the if statement */
      if (cmesh_id < 97 || cmesh_id > 256) {
        t8_test_emelemts_remove(cmesh_id);
      }
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  srand(time(0));
  //test_cmesh_emelemts_remove_all ();
  t8_test_emelemts_remove(0);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
