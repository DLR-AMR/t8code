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

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>

int
main (int argc, char **argv)
{
  /* The prefix for our output files. */
  const char prefix[BUFSIZ] = "t8_cmesh_hypercube_pad";
  t8_locidx_t local_num_trees;
  t8_gloidx_t global_num_trees;

  const double boundary_coords[24] = { 1, 0, 0, 4, 0, 0, 0, 6, 0, 5, 5, 0, -1, -2, 8, 9, 0, 10, 0, 8, 9, 10, 10, 10 };

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Add hypercube with given element class. */
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, boundary_coords, 3, 3, 3);

  /* Compute local and global number of trees. */
  local_num_trees = t8_cmesh_get_num_local_trees (cmesh);
  global_num_trees = t8_cmesh_get_num_trees (cmesh);
  t8_global_productionf (" [step1] Created coarse mesh.\n");
  t8_global_productionf (" [step1] Local number of trees:\t%i\n", local_num_trees);
  t8_global_productionf (" [step1] Global number of trees:\t%li\n", global_num_trees);

  t8_cmesh_vtk_write_file (cmesh, prefix, 1.0);

  t8_cmesh_destroy (&cmesh);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
