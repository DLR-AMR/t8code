/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

#include <test/t8_gtest_macros.hxx>

/**
 * Package id for the testsuite. Used for attributes.
 */
static int testsuite_package_id = -1;

void
t8_testsuite_init (int *argc, char ***argv, int log_threshold)
{
  /* Initialize mpi */
  const int mpiret = sc_MPI_Init (argc, argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc and t8code */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (log_threshold);

  /* Register a package id for the t8code testsuite */
  testsuite_package_id = sc_package_register (
    NULL, SC_LP_DEFAULT, "t8code_testsuite", "t8code testsuite package. Used for testing of external user attributes.");
}

void
t8_testsuite_finalize ()
{
  /* Finalize SC */
  sc_finalize ();

  /* Finalize and check mpi */
  const int mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
}

int
t8_testsuite_get_package_id ()
{
  return testsuite_package_id;
}
