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

#include <sc_refcount.h>
#include <t8_default.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <sc_shmem.h>

/* Create a partition that concentrates everything at a given proc */
static t8_gloidx_t *
t8_partition_offset (int proc, sc_MPI_Comm comm, t8_gloidx_t num_trees)
{
  int               mpirank, mpiret, mpisize, iproc;
  t8_gloidx_t      *offsets;
#ifdef T8_ENABLE_DEBUG
  char              out[BUFSIZ] = "";
#endif

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  offsets = SC_SHMEM_ALLOC (t8_gloidx_t, mpisize + 1,
                            sc_MPI_COMM_WORLD);
  offsets[0] = 0;
  for (iproc = 1; iproc <= mpisize; iproc++) {
    if (iproc == proc + 1) {
      offsets[iproc] = num_trees;
    }
    else {
      offsets[iproc] = offsets[iproc - 1];
    }
#ifdef T8_ENABLE_DEBUG
    snprintf (out + strlen(out), BUFSIZ - strlen(out), "%li,",offsets[iproc]);
#endif
  }
#ifdef T8_ENABLE_DEBUG
  t8_debugf ("Partition with offsets:0,%s\n", out);
#endif
  return offsets;
}

static void
t8_partition ()
{
  t8_cmesh_t        cmesh, cmesh_part, cmesh_part2;
  char              file[BUFSIZ];
  p4est_connectivity_t *conn;
  int               mpirank, mpiret, i, mpisize;


  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  conn = p4est_connectivity_new_brick (2,2,0,0);
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, 0, 1);
  p4est_connectivity_destroy (conn);
  snprintf (file, BUFSIZ,"t8_brick_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh, file, 1.);

  t8_cmesh_init (&cmesh_part);
  t8_cmesh_set_partition_from (cmesh_part, cmesh, -1,
                               t8_partition_offset (1, sc_MPI_COMM_WORLD,
                               t8_cmesh_get_num_trees (cmesh)));
  t8_cmesh_commit (cmesh_part);
  t8_cmesh_init (&cmesh_part2);
  t8_cmesh_set_partition_from (cmesh_part2, cmesh_part, -1,
                               t8_partition_offset (0, sc_MPI_COMM_WORLD,
                               t8_cmesh_get_num_trees (cmesh)));
  t8_cmesh_commit (cmesh_part2);

  snprintf (file, BUFSIZ,"t8_brick_partition_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh_part, file, 1.0);
  snprintf (file, BUFSIZ,"t8_brick_partition2_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh_part2, file, 1.0);
  t8_cmesh_unref (&cmesh_part2);
}



int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 level;
  int                 eclass;
  int                 dim;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_partition ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
