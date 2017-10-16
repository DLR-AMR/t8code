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

static void
test_cmesh_copy (sc_MPI_Comm comm)
{
  int                 eci, retval;
  t8_cmesh_t          cmesh_original, cmesh_copy;

  for (eci = T8_ECLASS_ZERO; eci < T8_ECLASS_COUNT; ++eci) {
    t8_global_productionf ("Testing eclass %s.\n", t8_eclass_to_string[eci]);

    /* Create new hypercube cmesh */
    cmesh_original = t8_cmesh_new_hypercube (eci, comm, 0, 0, 0);
    test_cmesh_committed (cmesh_original);
    /* Set up the cmesh copy */
    t8_cmesh_init (&cmesh_copy);
    /* We need the original cmesh later, so we ref it */
    t8_cmesh_ref (cmesh_original);
    t8_cmesh_set_derive (cmesh_copy, cmesh_original);
    /* Commit and check commit */
    t8_cmesh_commit (cmesh_copy, comm);
    test_cmesh_committed (cmesh_copy);
    /* Check for equality */
    retval = t8_cmesh_is_equal (cmesh_copy, cmesh_original);
    SC_CHECK_ABORT (retval == 1, "Cmesh copy failed.");
    /* Clean-up */
    t8_cmesh_destroy (&cmesh_copy);
    t8_cmesh_destroy (&cmesh_original);
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
  test_cmesh_copy (comm);
  t8_global_productionf ("Done testing cmesh copy.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
