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
#include <t8_cmesh_readmshfile.h>

/* TODO: rename this file to t8_something */

static void
t8_cmesh_save_cmesh (char *mshfile)
{
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  int                 ret, mpirank, mpiret;

  if (mshfile == NULL) {
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_WORLD,
                                    0, 0, 1);
  }
  else {
    t8_cmesh_t          cmesh_partition;
    cmesh = t8_cmesh_from_msh_file (mshfile, 1, sc_MPI_COMM_WORLD, 3, 0);
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_partition, 1);
    t8_cmesh_commit (cmesh_partition, sc_MPI_COMM_WORLD);
    t8_cmesh_destroy (&cmesh);
    cmesh = cmesh_partition;
  }
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  snprintf (filename, BUFSIZ, "cmesh_saved_%04d.cmsh", mpirank);
  ret = t8_cmesh_save (cmesh, filename);
  if (ret == 0) {
    t8_errorf ("Error when writing to file\n");
  }
  else {
    t8_debugf ("Saved cmesh to %s\n", filename);
  }
  t8_cmesh_destroy (&cmesh);
}

static void
t8_cmesh_load_cmesh ()
{
  t8_cmesh_t          cmesh;
  char                filename[BUFSIZ];
  int                 mpirank, mpiret;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  snprintf (filename, BUFSIZ, "cmesh_saved_%04d.cmsh", mpirank);
  cmesh = t8_cmesh_load (filename, sc_MPI_COMM_WORLD);
  if (cmesh != NULL) {
    t8_debugf ("Successfully loaded cmesh from %s\n", filename);
    t8_cmesh_vtk_write_file (cmesh, "cmesh_loaded", 1.0);
    t8_cmesh_destroy (&cmesh);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  char               *meshfile;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  if (argc == 2) {
    meshfile = argv[1];
  }
  else {
    meshfile = NULL;
  }
  t8_cmesh_save_cmesh (meshfile);
  t8_cmesh_load_cmesh ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
