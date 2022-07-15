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
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_testcases.h>

/* Test if a cmesh is committed properly and perform the
 * face consistency check. */

static void
t8_test_cmesh_committed (t8_cmesh_t cmesh)
{
  int                 retval;

  retval = t8_cmesh_is_committed (cmesh);
  SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
  retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
  SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
}

/** The function test_cmesh_copy (int cmesh_id,sc_MPI_Comm comm) runs the cmesh_copy test for one given cmesh,
 * that we get through its id by caling t8_test_create_cmesh (cmesh_id). 
 * \param [in] cmesh_id The cmesh_id which is used to create a unique cmesh with t8_test_create_cmesh.
 * \param [in] comm The communicator used to commit the cmesh_copy.
 */
static void
test_cmesh_copy (int cmesh_id, sc_MPI_Comm comm)
{
  int                 retval;
  t8_cmesh_t          cmesh_original, cmesh_copy;
  /* Create new cmesh */
  cmesh_original = t8_test_create_cmesh (cmesh_id);
  t8_test_cmesh_committed (cmesh_original);
  /* Set up the cmesh copy */
  t8_cmesh_init (&cmesh_copy);
  /* We need the original cmesh later, so we ref it */
  t8_cmesh_ref (cmesh_original);
  t8_cmesh_set_derive (cmesh_copy, cmesh_original);
  /* Commit and check commit */
  t8_cmesh_commit (cmesh_copy, comm);
  t8_test_cmesh_committed (cmesh_copy);
  /* Check for equality */
  retval = t8_cmesh_is_equal (cmesh_copy, cmesh_original);
  SC_CHECK_ABORT (retval == 1, "Cmesh copy failed.");
  /* Clean-up */
  t8_cmesh_destroy (&cmesh_copy);
  t8_cmesh_destroy (&cmesh_original);
}

/** The function test_cmesh_copy_all(sc_MPI_Comm comm) runs the cmesh_copy test for all cmeshes we want to test.
 * We run over all testcases using t8_get_all_testcases() to know how many to check. 
 * \param [in] comm The communicator used in test_cmesh_copy to commit the cmesh copy.
 */
static void
test_cmesh_copy_all (sc_MPI_Comm comm)
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases ();
       cmesh_id++) {
    test_cmesh_copy (cmesh_id, comm);
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

  t8_global_productionf ("Testing cmesh copy.\n");
  test_cmesh_copy_all (comm);

  t8_global_productionf ("Done testing cmesh copy.\n");
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
