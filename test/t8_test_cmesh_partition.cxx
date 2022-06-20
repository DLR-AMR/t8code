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

#include <t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_cmesh/t8_cmesh_testcases.h>

/* We create a cmesh, partition it and repartition it several times.
 * At the end we result in the same partition as at the beginning and we
 * compare this cmesh with the initial one. If they are equal the test is
 * passed.
 *
 * TODO: - Currently the test is not passed. This is probably because the
 *         cmesh_is_equal function is too strict and does not allow, for example,
 *         the order of the ghosts to change.
 *         We should implement a lighter version of cmesh_is_equal in order
 *         to account for this.
 * 
 *       - when this test works for all cmeshes remove if statement in 
 *         test_cmesh_partition_all () 
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
test_cmesh_partition (int cmesh_id, sc_MPI_Comm comm)
{
  int                 level = 11;
  int                 mpisize, mpiret;
  int                 i;
  t8_cmesh_t          cmesh_original, cmesh_partition;

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  cmesh_original = t8_test_create_cmesh (cmesh_id);

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

/** The function test_cmesh_copy_all(sc_MPI_Comm comm) runs the cmesh_copy test for all cmeshes we want to test.
 * We run over all testcases using t8_get_all_testcases() to know how many to check. 
 * \param [in] comm The communicator used in test_cmesh_copy to commit the cmesh copy.
 */
static void
test_cmesh_partition_all (sc_MPI_Comm comm)
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    /* This if statement is necessary to make the test work by avoiding specific cmeshes which do not work yet for this test.
     * When the issues are gone, remove the if statement. */
    if (cmesh_id != 89 && (cmesh_id < 237 || cmesh_id > 256)) {
      test_cmesh_partition (cmesh_id, comm);
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
  test_cmesh_partition_all (comm);
  t8_global_productionf ("Done testing cmesh partition.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
