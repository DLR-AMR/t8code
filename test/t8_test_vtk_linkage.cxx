/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* In this test we create a vtk unstructured Grid object and check whether
 * vtk version t8code was configured to link with is the version that is actually
 * linked.
 * The purpose of this test is to check whether t8code successfully links
 * against vtk.
 * If t8code was not configured with --with-vtk then this test
 * does nothing and is always passed.
 */

#include <t8.h>
#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>
#include <vtkNew.h>
#endif

/* Test correct macro dependencies.
 * Will throw a compile time error if T8_VTK_VERSION_USED
 * is defined but T8_WITH_VTK is not. */
#ifndef T8_WITH_VTK
#ifdef T8_VTK_VERSION_USED
#error Configuration error: T8_VTK_VERSION_USED is defined despite \
 T8_WITH_VTK not being defined.
#endif
#endif

/* Check whether T8_VTK_VERSION_USED equals VTK_MAJOR_VERSION.VTK_MINOR_VERSION */
static void
t8_test_vtk_version_number ()
{
#if T8_WITH_VTK
  char                vtk_version[BUFSIZ];
  snprintf (vtk_version, BUFSIZ, "%i.%i", VTK_MAJOR_VERSION,
            VTK_MINOR_VERSION);
  if (strcmp (T8_VTK_VERSION_USED, vtk_version)) {
    SC_ABORTF
      ("linked vtk version (%s) does not equal the version t8code was configured"
       " with (%s)\n", vtk_version, T8_VTK_VERSION_USED);
  }
  else {
    t8_global_productionf ("Using vtk version %s.\n", vtk_version);
  }
#endif
}

/* Check whether we can successfully execute VTK code */
static void
t8_test_vtk_linkage ()
{
#if T8_WITH_VTK
  vtkNew < vtkUnstructuredGrid > unstructuredGrid;

  t8_global_productionf
    ("Successfully created VTK unstructuredGrid object.\n");
#else
  t8_global_productionf
    ("This version of t8code is not compiled with vtk support.\n");
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  /* Check for correct version number */
  t8_test_vtk_version_number ();
  /* Check vtk linkage */
  t8_test_vtk_linkage ();

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
