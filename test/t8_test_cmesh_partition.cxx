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

#include <t8_cmesh.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"

/* We create a cmesh, partition it and repartition it several times.
 * At the end we result in the same partition as at the beginning and we
 * compare this cmesh with the initial one. If they are equal the test is
 * passed.
 *
 * TODO: Currently the test is not passed. This is probably because the
 *       cmesh_is_equal function is too strict and does not allow, for example,
 *       the order of the ghosts to change.
 *       We should implement a lighter version of cmesh_is_equal in order
 *       to account for this.
 */

/* Test if a cmesh is committed properly and perform the
 * face consistency check. */
static void
test_cmesh_committed (t8_cmesh_t cmesh)
{
  int                 retval;

  retval = t8_cmesh_is_committed (cmesh);
  SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
  retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
  SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
}

/* For each process re-partition a partitioned cmesh to be concentrated on
 * this process and check for equality with the non-partitioned cmesh.
 */
static void
test_cmesh_partition_concentrate (t8_cmesh_t cmesh_partition_uniform,
                                  sc_MPI_Comm comm, int level)
{
  int                 mpisize, mpirank, mpiret, irank;
  int                 retval, i;
  t8_cmesh_t          cmesh_partition_new1, cmesh_partition_new2;
  t8_shmem_array_t    offset_concentrate;

  /* Get the number of processes and rank of current process */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Since we want to repartition the cmesh_partition in each step,
   * we need to ref it. This ensure that we can still work with it after
   * another cmesh is derived from it. */
  t8_cmesh_ref (cmesh_partition_uniform);

  cmesh_partition_new1 = cmesh_partition_uniform;

  /* We repartition the cmesh to be concentrated on each rank once */
  for (irank = 0; irank < mpisize; irank++) {
    t8_cmesh_init (&cmesh_partition_new2);
    t8_cmesh_set_derive (cmesh_partition_new2, cmesh_partition_new1);
    /* Create an offset array where each tree resides on irank */
    offset_concentrate = t8_cmesh_offset_concentrate (irank, comm,
                                                      t8_cmesh_get_num_trees
                                                      (cmesh_partition_uniform));
    /* Set the new cmesh to be partitioned according to that offset */
    t8_cmesh_set_partition_offsets (cmesh_partition_new2, offset_concentrate);
    /* Commit the cmesh and test if successful */
    t8_cmesh_commit (cmesh_partition_new2, comm);
    test_cmesh_committed (cmesh_partition_new2);

    /* Switch the rolls of the cmeshes */
    cmesh_partition_new1 = cmesh_partition_new2;
    cmesh_partition_new2 = NULL;
  }
  /* We partition the resulting cmesh according to a uniform level refinement.
   * This cmesh should now be equal to the initial cmesh. */
  for (i = 0; i < 2; i++) {
    t8_cmesh_init (&cmesh_partition_new2);
    t8_cmesh_set_derive (cmesh_partition_new2, cmesh_partition_new1);
    t8_cmesh_set_partition_uniform (cmesh_partition_new2, level,
                                    t8_scheme_new_default_cxx ());
    t8_cmesh_commit (cmesh_partition_new2, comm);
    cmesh_partition_new1 = cmesh_partition_new2;
  }
  retval = t8_cmesh_is_equal (cmesh_partition_new2, cmesh_partition_uniform);
  SC_CHECK_ABORT (retval == 1, "Cmesh equality check failed.");

  /* clean-up */
  t8_cmesh_destroy (&cmesh_partition_new2);
}

static void
test_cmesh_partition (sc_MPI_Comm comm)
{
  int                 eci, level, maxlevel, minlevel;
  int                 mpisize, mpiret;
  int                 i;
  t8_cmesh_t          cmesh_original, cmesh_partition;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  for (eci = T8_ECLASS_VERTEX; eci < T8_ECLASS_COUNT; ++eci) {
    t8_global_productionf ("\n\nTesting eclass %s.\n",
                           t8_eclass_to_string[eci]);

    if (eci != T8_ECLASS_PYRAMID) {
      /* TODO: cmesh_partition does not work with pyramids yet.
       *       as soon as it does, activate the test for pyramids.
       */

      /* We do not support empty forest processes for cmesh uniform partition
       * yet, thus we must ensure that the number of forest elements in a
       * uniform refinement would be equal to or exceed the number of processes.
       * We compute here the number of forest elements per level and set as
       * minimum refinement level the first level where the condition is
       * fulfilled. */
      /* For eclass vertex the forest always has 1 element when using the hypercube mesh
       * and thus we use a disjoint copy of at least mpisize many coarse cells. */

      minlevel = 0;
      maxlevel = minlevel + 8;
      for (level = minlevel; level < maxlevel; level++) {
        t8_global_productionf ("\n\tTesting refinement level %i\n", level);
        if (eci != T8_ECLASS_VERTEX && level != maxlevel - 1) {
          /* Take the hypercube as coarse mesh */
          cmesh_original =
            t8_cmesh_new_hypercube ((t8_eclass_t) eci, comm, 0, 0, 0);
        }
        else {
          /* If the eclass is vertex we choose a mesh consisting of disjoint
           * vertices.
           * We also choose this mesh for all eclasses in the last case.
           * We do this to test different mesh sizes. */
          double              num_trees = 11;   /* prime number */
          cmesh_original =
            t8_cmesh_new_bigmesh ((t8_eclass_t) eci, num_trees, comm);
        }
        test_cmesh_committed (cmesh_original);
        for (i = 0; i < 2; i++) {
          /* Set up the partitioned cmesh */
          t8_cmesh_init (&cmesh_partition);
          t8_cmesh_set_derive (cmesh_partition, cmesh_original);
          /* Uniform partition according to level */
          t8_cmesh_set_partition_uniform (cmesh_partition, level,
                                          t8_scheme_new_default_cxx ());
          t8_cmesh_commit (cmesh_partition, comm);
          test_cmesh_committed (cmesh_partition);
          cmesh_original = cmesh_partition;
        }

        /* Perform the concentrate test */
        test_cmesh_partition_concentrate (cmesh_partition, comm, level);
        /* Clean-up */
        t8_cmesh_destroy (&cmesh_partition);
      }
    }
    else {
      t8_global_productionf ("Skipping test for %s.\n",
                             t8_eclass_to_string[eci]);
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  comm = sc_MPI_COMM_WORLD;
  sc_init (comm, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf ("Testing cmesh partition.\n");
  test_cmesh_partition (comm);
  t8_global_productionf ("Done testing cmesh partition.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
