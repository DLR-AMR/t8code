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
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_default_cxx.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>

/* In this test, we adapt, balance and partition a uniform forest.
 * We do this in two ways:
 * 1st  All operations are performed in one single call to t8_forest_commit
 * 2nd  Each intermediate step is performed in a seperate commit
 *
 * After these two forests are created, we check for equality.
 */

/* Adapt a forest such that always the second child of a
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_balance (t8_forest_t forest, t8_forest_t forest_from,
                       t8_locidx_t which_tree, t8_eclass_scheme_c * ts,
                       int num_elements, t8_element_t * elements[])
{
  int                 level;
  int                 maxlevel, child_id;
  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));
  level = ts->t8_element_level (elements[0]);

  /* we set a maximum refinement level as forest user data */
  maxlevel = *(int *) t8_forest_get_user_data (forest);
  if (level >= maxlevel) {
    /* Do not refine after the maxlevel */
    return 0;
  }
  child_id = ts->t8_element_child_id (elements[0]);
  if (child_id == 2) {
    return 1;
  }
  return 0;
}

/* Depending on an integer i create a different cmesh.
 * i = 0: cmesh_new_class
 * i = 1: cmesh_new_hypercube
 * i = 2: cmesh_new_bigmesh (100 trees)
 * else:  cmesh_new_class
 */
static              t8_cmesh_t
t8_test_create_cmesh (int i, t8_eclass_t eclass, sc_MPI_Comm comm)
{
  switch (i) {
  case 0:
    return t8_cmesh_new_from_class (eclass, comm);
  case 1:
    return t8_cmesh_new_hypercube (eclass, comm, 0, 0);
  case 2:
    return t8_cmesh_new_bigmesh (eclass, 100, comm);
  default:
    return t8_cmesh_new_from_class (eclass, comm);
  }
}

/* adapt, balance and partition a given forest in one step */
static              t8_forest_t
t8_test_forest_commit_abp (t8_forest_t forest, int maxlevel)
{
  t8_forest_t         forest_ada_bal_par;

  /* Adapt, balance and partition the uniform forest */
  t8_forest_init (&forest_ada_bal_par);
  /* Set user data for adapt */
  t8_forest_set_user_data (forest_ada_bal_par, &maxlevel);
  t8_forest_set_adapt (forest_ada_bal_par, forest, t8_test_adapt_balance,
                       NULL, 1);
  t8_forest_set_balance (forest_ada_bal_par, NULL, 0);
  t8_forest_set_partition (forest_ada_bal_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_par);

  t8_debugf ("[H] 1 step forest %i local %li global elements.\n",
             t8_forest_get_num_element (forest_ada_bal_par),
             t8_forest_get_global_num_elements (forest_ada_bal_par));

  return forest_ada_bal_par;
}

/* adapt, balance and partition a given forest in 3 steps */
static              t8_forest_t
t8_test_forest_commit_abp_3step (t8_forest_t forest, int maxlevel)
{
  t8_forest_t         forest_adapt, forest_balance, forest_partition;

  t8_forest_init (&forest_adapt);
  t8_forest_init (&forest_balance);
  t8_forest_init (&forest_partition);

  /* adapt the forest */
  t8_forest_set_user_data (forest_adapt, &maxlevel);
  t8_forest_set_adapt (forest_adapt, forest, t8_test_adapt_balance, NULL, 1);
  t8_forest_commit (forest_adapt);

  /* balance the forest */
  t8_forest_set_balance (forest_balance, forest_adapt, 0);
  t8_forest_commit (forest_balance);

  /* partrition the forest */
  t8_forest_set_partition (forest_partition, forest_balance, 0);
  t8_forest_commit (forest_partition);
  t8_debugf ("[H] 3 step forest %i local %li global elements.\n",
             t8_forest_get_num_element (forest_partition),
             t8_forest_get_global_num_elements (forest_partition));

  return forest_partition;
}

static void
t8_test_forest_commit ()
{
  int                 ctype, level, min_level, maxlevel;
  int                 eclass;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest, forest_ada_bal_part, forest_abp_3part;
  t8_scheme_cxx_t    *scheme;

  for (eclass = T8_ECLASS_QUAD; eclass < T8_ECLASS_PRISM; eclass++) {
    /* TODO: Activate the other eclass as soon as they support ghosts */
    //for (ctype = 0; ctype < 3; ctype++) {
    {
      scheme = t8_scheme_new_default_cxx ();
      ctype = 2;
      /* Construct a cmesh */
      cmesh =
        t8_test_create_cmesh (ctype, (t8_eclass_t) eclass, sc_MPI_COMM_WORLD);
      min_level = t8_forest_min_nonempty_level (cmesh, scheme);
      t8_global_productionf
        ("Testing forest commit with eclass %s, start level %i\n",
         t8_eclass_to_string[eclass], min_level);
      for (level = min_level; level < min_level + 3; level++) {
        maxlevel = level + 3;
        /* ref the scheme since we reuse it */
        t8_scheme_cxx_ref (scheme);
        /* ref the cmesh since we reuse it */
        t8_cmesh_ref (cmesh);
        /* Create a uniformly refined forest */
        forest = t8_forest_new_uniform (cmesh, scheme, level, 1,
                                        sc_MPI_COMM_WORLD);
        /* We need to use forest twice, so we ref it */
        t8_forest_ref (forest);
        /* Adapt, balance and partition the forest */
        forest_ada_bal_part = t8_test_forest_commit_abp (forest, maxlevel);
        /* Adapt, balance and partition the forest using three seperate steps */
        forest_abp_3part = t8_test_forest_commit_abp_3step (forest, maxlevel);
        if (ctype != 2) {
          t8_forest_write_vtk (forest_ada_bal_part, "test_1step");
          t8_forest_write_vtk (forest_abp_3part, "test_3step");
        }
        SC_CHECK_ABORT (t8_forest_is_equal
                        (forest_abp_3part, forest_ada_bal_part),
                        "The forests are not equal");

        t8_forest_unref (&forest_ada_bal_part);
        t8_forest_unref (&forest_abp_3part);
      }
      t8_cmesh_destroy (&cmesh);
      t8_scheme_cxx_unref (&scheme);
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

  t8_test_forest_commit ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
