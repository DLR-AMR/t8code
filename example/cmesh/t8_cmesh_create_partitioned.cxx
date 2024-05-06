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

#include <t8_cmesh_vtk_writer.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>

void
t8_cmesh_create_partitioned (sc_MPI_Comm comm)
{
  int mpisize, mpirank;
  t8_cmesh_t cmesh;

  t8_cmesh_init (&cmesh);

  /* Get the number of processes */
  sc_MPI_Comm_size (comm, &mpisize);
  /* Get this processes rank */
  sc_MPI_Comm_rank (comm, &mpirank);

  /* Set the classes of the trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);

  /* Set the face-to-face connections of the owned trees */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0); /* tree 0 - tree 1 */

  t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
  t8_cmesh_set_partition_uniform (cmesh, 0, scheme);

  /* Create the cmesh */
  t8_cmesh_commit (cmesh, comm);

  /* Clean-up */
  t8_cmesh_destroy (&cmesh);
}

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_cmesh_create_partitioned (sc_MPI_COMM_WORLD);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
