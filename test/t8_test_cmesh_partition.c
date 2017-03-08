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
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"

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
test_cmesh_partition_concentrate (t8_cmesh_t cmesh_original,
                                  t8_cmesh_t cmesh_partition,
                                  sc_MPI_Comm comm)
{
  int                 mpisize, mpirank, mpiret, irank;
  int                 retval;
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
  t8_cmesh_ref (cmesh_partition);

  t8_cmesh_init (&cmesh_partition_new1);
  t8_cmesh_set_derive (cmesh_partition_new1, cmesh_partition);
  t8_cmesh_commit (cmesh_partition_new1, comm);
  t8_cmesh_unref (&cmesh_partition_new1);

#if 0
  for (irank = 0; irank < mpisize; irank++) {
    t8_cmesh_init (&cmesh_concentrate);
    t8_cmesh_set_derive (cmesh_concentrate, cmesh_partition);
    /* Create an offset array where each tree resides on irank */
    offset_concentrate = t8_cmesh_offset_concentrate (irank, comm,
                                                      t8_cmesh_get_num_trees
                                                      (cmesh_partition));
    /* Set the new cmesh to be partitioned according to that offset */
    t8_cmesh_set_partition_offsets (cmesh_concentrate, offset_concentrate);
    /* Commit the cmesh and test if successful */
    t8_cmesh_commit (cmesh_concentrate, comm);
    test_cmesh_committed (cmesh_concentrate);
    if (mpirank == irank) {
      /* If all trees are concentrated on the current process, check whether
       * this cmesh is equal to the original cmesh */
      retval = t8_cmesh_is_equal (cmesh_concentrate, cmesh_original);
      SC_CHECK_ABORT (retval == 1, "Cmesh equality check failed.");
    }
    t8_cmesh_unref (&cmesh_concentrate);
    t8_cmesh_unref (&cmesh_partition);
  }
#endif
}

static void
test_cmesh_partition (sc_MPI_Comm comm)
{
  int                 eci, level;
  t8_cmesh_t          cmesh_original, cmesh_partition;

  for (eci = T8_ECLASS_ZERO; eci < T8_ECLASS_LINE; ++eci) {
    t8_global_productionf ("Testing eclass %s.\n", t8_eclass_to_string[eci]);
    if (eci != T8_ECLASS_PYRAMID) {
      /* TODO: cmesh_partition does not work with pyramids yet.
       *       as soon as it does, activate the test for pyramids.
       */
      for (level = 0; level < 4; level++) {
        cmesh_original = t8_cmesh_new_hypercube (eci, comm, 0, 0);
        test_cmesh_committed (cmesh_original);
        /* Set up the partitioned cmesh */
        t8_cmesh_init (&cmesh_partition);
        t8_debugf ("  original: %i\n", cmesh_original->rc.refcount);
        t8_cmesh_set_derive (cmesh_partition, cmesh_original);
        /* Uniform partition according to level */
        t8_cmesh_set_partition_uniform (cmesh_partition, level);
        t8_cmesh_commit (cmesh_partition, comm);
        t8_debugf ("  original: %i\n", cmesh_original->rc.refcount);
        test_cmesh_committed (cmesh_partition);
        /* Perform the concentrate test */
#if 0
        /* TODO: fix reference counting */
        test_cmesh_partition_concentrate (cmesh_original, cmesh_partition,
                                          comm);
#endif
        /* Clean-up */
        t8_cmesh_unref (&cmesh_partition);
        t8_debugf ("  original: %i\n", cmesh_original->rc.refcount);
        t8_cmesh_unref (&cmesh_original);
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
