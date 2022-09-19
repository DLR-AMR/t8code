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
#include <t8_cmesh/t8_cmesh_examples.h>

static void
test_hypercube (sc_MPI_Comm mpic)
{
  int                 eci;
  int                 retval;
  int                 partition, bcast;
  t8_cmesh_t          cmesh;

  for (eci = T8_ECLASS_ZERO; eci < T8_ECLASS_COUNT; ++eci) {
    for (partition = 0; partition < 2; partition++) {
      for (bcast = 0; bcast < 2; bcast++) {
        cmesh = t8_cmesh_new_hypercube
          ((t8_eclass_t) eci, mpic, bcast, partition, 0);
        retval = t8_cmesh_is_committed (cmesh);
        SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
        retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
        SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
        t8_cmesh_destroy (&cmesh);
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

  test_hypercube (mpic);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
