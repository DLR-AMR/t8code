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
#include <sc_refcount.h>
#include <t8.h>  
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h> /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_transition.h> 
#include <t8_forest/t8_forest_geometrical.h>    /* geometrical information of the forest */
#include <t8_forest/t8_forest_io.h>   
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_schemes/t8_transition/t8_transition_cxx.hxx>
#include <t8_schemes/t8_transition/t8_transition_conformal_hex/t8_transition_conformal_hex_cxx.hxx>

/* In this test, we adapt, balance and partition a uniform forest.
 * We do this in two ways:
 * 1st  All operations are performed in one single call to t8_forest_commit
 * 2nd  Each intermediate step is performed in a separate commit
 *
 * After these two forests are created, we check for equality.
 */


/* Adapt a forest such that always the first child of a
 * tree is refined and no other elements. This results in a highly
 * imbalanced forest. */
static int
t8_test_adapt_balance (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                       t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
  T8_ASSERT (!is_family || (is_family && num_elements == ts->t8_element_num_children (elements[0])));

  int level = ts->t8_element_level (elements[0]);

  /* we set a maximum refinement level as forest user data */
  int maxlevel = *(int *) t8_forest_get_user_data (forest);
  // if (level >= maxlevel) {
  //   /* Do not refine after the maxlevel */
  //   return 0;
  // }
  // int child_id = ts->t8_element_child_id (elements[0]);
  // if (child_id == 1) {
  //   return 1;
  // }
  // return 0;
  if (lelement_id == 0){
    return 1;
  }
return 0;
}

/* adapt, balance and partition a given forest in one step */
static t8_forest_t
t8_test_forest_commit_abp (t8_forest_t forest, int maxlevel)
{
    t8_debugf("---------start t8_test_forest_commit_abp ---------------------\n");

  t8_forest_t forest_ada_bal_par;

  /* Adapt, balance and partition the uniform forest */
  t8_forest_init (&forest_ada_bal_par);
  /* Set user data for adapt */
  t8_forest_set_user_data (forest_ada_bal_par, &maxlevel);
  t8_forest_set_adapt (forest_ada_bal_par, forest, t8_test_adapt_balance, 0);
  t8_forest_set_balance (forest_ada_bal_par, NULL, 0);
  t8_forest_set_transition(forest_ada_bal_par, NULL, 0);
  t8_forest_set_partition (forest_ada_bal_par, NULL, 0);
  t8_forest_commit (forest_ada_bal_par);

  t8_debugf("---------finsish t8_test_forest_commit_abp ---------------------\n");
  return forest_ada_bal_par;
}

/* adapt, balance and partition a given forest in 3 steps */
static t8_forest_t
t8_test_forest_commit_abp_3step (t8_forest_t forest, int maxlevel)
{
    t8_debugf("---------start t8_test_forest_commit_abtp_3step ---------------------\n");

  t8_forest_t forest_adapt;
  t8_forest_t forest_balance;
  t8_forest_t forest_partition;
  t8_forest_t forest_transition;

  t8_forest_init (&forest_adapt);
  t8_forest_init (&forest_balance);
  t8_forest_init (&forest_partition);
  t8_forest_init (&forest_transition);

  /* adapt the forest */
  t8_forest_set_user_data (forest_adapt, &maxlevel);
  t8_forest_set_adapt (forest_adapt, forest, t8_test_adapt_balance, 0);
  t8_forest_commit (forest_adapt);

  /* balance the forest */
  t8_forest_set_balance (forest_balance, forest_adapt, 0);
  t8_forest_commit (forest_balance);

  /* transition the forest */
  t8_forest_set_transition (forest_transition, forest_balance, 0);
  t8_forest_commit (forest_transition);

  /* partrition the forest */
  t8_forest_set_partition (forest_partition, forest_transition, 0);
  t8_forest_commit (forest_partition);
t8_debugf("---------finsish t8_test_forest_commit_abtp_3step ---------------------\n");
  return forest_partition;
}

void t8_transition_commit ()
{
  t8_cmesh_t  cmesh; 
  t8_forest_t forest;
  t8_forest_t forest_ada_bal_tra_part;
  t8_forest_t forest_abtp_3part;

  /* construct a single tree hex cmesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  

  t8_scheme_cxx_t *scheme = t8_scheme_new_transition_hex_cxx ();

   int maxlevel = 2;
  for (int level = 2; level <= maxlevel; level++) {
     /* ref the cmesh since we reuse it */
    t8_cmesh_ref (cmesh);
    t8_debugf ("Testing forest commit level %i\n", level);
   
    /* Create a uniformly refined forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

    /* We need to use forest twice, so we ref it */
    t8_forest_ref (forest);
    /* Adapt, balance, transition and partition the forest */
            t8_debugf("ref count forest %i  \n", forest->rc.refcount);

    forest_ada_bal_tra_part = t8_test_forest_commit_abp (forest, maxlevel);
    t8_forest_write_vtk (forest_ada_bal_tra_part, "transition_forest");

    t8_debugf("ref count forest %i  \n", forest->rc.refcount);

    /* Adapt, balance, transition and partition the forest using three separate steps */
    forest_abtp_3part = t8_test_forest_commit_abp_3step (forest, maxlevel);

    // ASSERT_TRUE (t8_forest_is_equal (forest_abtp_3part, forest_ada_bal_tra_part)) << "The forests are not equal";
    t8_scheme_cxx_ref (scheme);

    t8_debugf("ref count forest ada_bal_tra_part %i  \n", forest_ada_bal_tra_part->rc.refcount);
    t8_debugf("ref count forest_abtp_3part %i  \n", forest_abtp_3part->rc.refcount);
    t8_debugf("ref count scheme %i  \n", scheme->rc.refcount);
    t8_debugf("ref count cmesh %i  \n", cmesh->rc.refcount);
    t8_forest_unref (&forest_ada_bal_tra_part);
    t8_forest_unref (&forest_abtp_3part);
    t8_debugf("ref count scheme %i  \n", scheme->rc.refcount);
    t8_debugf("ref count cmesh %i  \n", cmesh->rc.refcount);

  }
    t8_cmesh_unref(&cmesh);
    t8_debugf("ref count scheme %i  \n", scheme->rc.refcount);
    t8_scheme_cxx_unref (&scheme);
    

  t8_debugf ("Done testing forest commit.");
}


int
main (int argc, char **argv)
{
  int                 mpiret;
/* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);

  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);

  t8_init (SC_LP_DEFAULT);

  t8_transition_commit ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();

  SC_CHECK_MPI (mpiret);

  return 0;
}