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

#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>

/* TODO: rename this file to t8_something */

void
t8_refine ()
{
  t8_cmesh_t          cmesh, cmesh_refine;

  cmesh =
    t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0, 0, 0);
  t8_cmesh_init (&cmesh_refine);
  t8_cmesh_set_derive (cmesh_refine, cmesh);
  t8_cmesh_set_refine (cmesh_refine, 1);
  t8_cmesh_commit (cmesh_refine, sc_MPI_COMM_WORLD);
  t8_cmesh_destroy (&cmesh_refine, sc_MPI_COMM_WORLD);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_refine ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
