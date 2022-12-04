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

#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "t8_cmesh/t8_cmesh_testcases.h"
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>

static int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t *elements[])
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
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t *elements[])
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
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t *elements[])
{
  int level_cut;
  int level = ts->t8_element_level (elements[0]);
  int level_max = ts->t8_element_maxlevel();
  
#if T8_ENABLE_LESS_TESTS
  level_cut = (int)(0.2*level_max);
#else
  level_cut = (int)(0.3*level_max);
#endif

  if (rand()%4 > 0 && level < level_cut) {
    return 1;
  }
  if (rand()%4 == 0) {
    return -2;
  }
  return 0;
}

static t8_forest_t
t8_adapt_forest (t8_forest_t forest_from,
                 t8_forest_adapt_t adapt_fn,
                 int do_partition,
                 int recursive,
                 int do_face_ghost)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  t8_forest_set_ghost (forest_new, do_face_ghost, T8_GHOST_FACES);
  if (do_partition) {
    t8_forest_set_partition (forest_new, NULL, 0);
  }
  t8_forest_commit (forest_new);

  return forest_new;
}

void
t8_test_elements_remove (int cmesh_id)
{
  int                 level;
  int                 level_min;
  int                 level_max;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();
  /* Construct a cmesh */
  cmesh = t8_test_create_cmesh (cmesh_id);
  /* Compute the first level, such that no process is empty */
  level_min = t8_forest_min_nonempty_level (cmesh, scheme);

#if T8_ENABLE_LESS_TESTS
  level_min = SC_MAX (level_min, 0);
  level_max = level_min + 1;
#else
  level_min = SC_MAX (level_min, 0);
  level_max = level_min + 2;
#endif

  for (level = level_min; level < level_max; level++) {
    t8_cmesh_ref (cmesh);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    for (int i = 0; i < level + 1; i++) {
        forest = t8_adapt_forest (forest, t8_adapt_callback_refine, 0, 0, 0);
        forest = t8_adapt_forest (forest, t8_adapt_callback_remove, 0, 0, 0);
    }
    forest = t8_adapt_forest (forest, t8_adapt_callback_refine, 0, 1, 0);
    for (int i = 0; i < 2*level_max; i++) {
      forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 0, 0, 0);
      forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 0, 0, 0);
      SC_CHECK_ABORT (t8_forest_no_overlap (forest),
          "The local forest has overlapping elements.");
      forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 0, 1, 0);
      SC_CHECK_ABORT (t8_forest_no_overlap (forest),
          "The local forest has overlapping elements.");
    }
    t8_scheme_cxx_ref (scheme);
    t8_forest_unref (&forest);
  }
  t8_scheme_cxx_unref (&scheme);
  t8_cmesh_destroy (&cmesh);
}

void
t8_test_cmesh_elements_remove_all ()
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases (); cmesh_id++) {
    t8_test_elements_remove (cmesh_id);
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
  t8_init (SC_LP_DEFAULT);

  unsigned int seed;
  seed = time(0);
  srand(seed);
  t8_global_productionf("Seed: %u \n", seed);

  t8_test_cmesh_elements_remove_all ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
