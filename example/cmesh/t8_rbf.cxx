/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_rbf.cxx
 *  This file implements an example for CAD-based mesh deformation.
 */
#include <t8.h>
#if T8_ENABLE_EIGEN
#include <Eigen/Dense>
#endif

int
main ([[maybe_unused]] int argc, [[maybe_unused]] char **argv)
{
#if T8_ENABLE_EIGEN
  t8_global_productionf ("--- STARTING EIGEN BASIC TEST --- \n");
  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);

  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);

  /* Initialize t8code with log level SC_LP_ESSENTIAL. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /** Little Test. */
  Eigen::Matrix3d A;
  A << 4.0, 1.2, 0.5, 1.2, 5.0, 2.1, 0.5, 2.1, 6.0;

  /** vector b. */
  Eigen::Vector3d b;
  b << 1.0, 2.0, 3.0;

  /* We will test a dense matrix example. Given that the matrix is symmetric and positive semidefinite we will test out the LDLT decomposition. */
  Eigen::LDLT<Eigen::Matrix3d> solver (A);

  if (solver.info () == Eigen::Success) {
    Eigen::Vector3d x = solver.solve (b);

    t8_global_productionf ("The solution is: x = [%f, %f, %f]\n", x (0), x (1), x (2));

    /* Test if A * x = b for a quick check. */
    Eigen::Vector3d check = A * x;
    t8_global_productionf ("(A*x): [%f, %f, %f] (should be the same as b)\n", check (0), check (1), check (2));
  }
  else {
    t8_global_productionf ("ERROR: the system could not be solved.\n");
  }

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

#else
  t8_global_productionf ("ERROR: This example requires Eigen support to be enabled in t8code.\n");
#endif

  return 0;
}
