
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

#include <time.h>
#include <sc_refcount.h>
#include <t8_default.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <sc_shmem.h>

static int
t8_offset_consistent (int mpisize, t8_gloidx_t * offset,
                      t8_gloidx_t num_trees)
{
  int                 i, ret = 1;
  t8_gloidx_t         temp;

  ret = offset[0] == 0;
  temp = 0;
  for (i = 1; i < mpisize && ret; i++) {
    if (offset[i] < 0) {
      ret &= (fabs (offset[i] + 1) >= temp);
      temp = fabs (offset[i] + 1);
    }
    else {
      ret &= (offset[i] >= temp);
      temp = offset[i];
    }
    ret &= (temp <= num_trees);
  }
  ret &= (offset[mpisize] == num_trees);
  return ret;
}

/* Create a partition that concentrates everything at a given proc */
static t8_gloidx_t *
t8_partition_offset (int proc, sc_MPI_Comm comm, t8_gloidx_t num_trees)
{
  int                 mpirank, mpiret, mpisize, iproc;
  t8_gloidx_t        *offsets;
#ifdef T8_ENABLE_DEBUG
  char                out[BUFSIZ] = "";
#endif

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  offsets = SC_SHMEM_ALLOC (t8_gloidx_t, mpisize + 1, sc_MPI_COMM_WORLD);
  offsets[0] = 0;
  for (iproc = 1; iproc <= mpisize; iproc++) {
    if (iproc == proc + 1) {
      offsets[iproc] = num_trees;
    }
    else {
      offsets[iproc] = offsets[iproc - 1];
    }
#ifdef T8_ENABLE_DEBUG
    snprintf (out + strlen (out), BUFSIZ - strlen (out), "%li,",
              offsets[iproc]);
#endif
  }
#ifdef T8_ENABLE_DEBUG
  t8_debugf ("Partition with offsets:0,%s\n", out);
#endif
  offsets[proc + 1] = -num_trees;
  return offsets;
}

/* Create a random partition */
/* if shared is nonzero than first trees can be shared */
static t8_gloidx_t *
t8_partition_offset_random (sc_MPI_Comm comm, t8_gloidx_t num_trees,
                            int shared)
{
  unsigned            seed;
  int                 iproc, mpisize, mpiret, random_number, mpirank;
  int                 first_shared;
  int                 i;
  t8_gloidx_t        *offsets, trees_so_far = 0;

  T8_ASSERT (shared == 0);      /* shared not yet implemented */

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  offsets = SC_SHMEM_ALLOC (t8_gloidx_t, mpisize + 1, sc_MPI_COMM_WORLD);

  do {
    seed = time (NULL);

    if (mpirank == 0) {
      t8_debugf ("Random number seed = %u\n", seed);
    }
    mpiret = sc_MPI_Bcast (&seed, 1, sc_MPI_INT, 0, sc_MPI_COMM_WORLD);
    SC_CHECK_MPI (mpiret);
    srand (seed);

    offsets[0] = 0;
    first_shared = 0;
    for (iproc = 1; iproc < mpisize; iproc++) {
      offsets[iproc] = 0;
      /* Create a random number between 0 and 200% of an ideal partition */
      random_number = rand () % (int)(num_trees * 2./mpisize);
      /* If we would excees the number of trees we cut the random number */
      if (offsets[iproc - 1] + random_number > num_trees) {
          random_number = num_trees - offsets[iproc - 1];
      }
      random_number += first_shared;
      first_shared = rand_r (&seed) % 2;

      offsets[iproc] = random_number + fabs (offsets[iproc - 1])
        - (offsets[iproc - 1] < 0);
      if (first_shared && offsets[iproc] != num_trees) {
        offsets[iproc] = -offsets[iproc] - 1;
      }
      trees_so_far += random_number;
    }
    offsets[mpisize] = num_trees;

    for (i = 0; i < mpisize; i++) {
      sc_MPI_Barrier (sc_MPI_COMM_WORLD);
      if (i == mpirank) {
        for (iproc = 0; iproc < mpisize + 1; iproc++) {
          t8_debugf ("[H] %li\n", offsets[iproc]);
        }
      }
    }
  }
  while (!t8_offset_consistent (mpisize, offsets, num_trees));
  return offsets;
}

static void
t8_random_partition ()
{
  t8_cmesh_t          cmesh, cmesh_part, cmesh_part2;
  char                file[BUFSIZ];
  p8est_connectivity_t *conn;
  int                 mpirank, mpiret, mpisize;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  conn = p8est_connectivity_new_brick (8, 8, 8, 0, 0, 0);

  cmesh = t8_cmesh_new_from_p8est (conn, sc_MPI_COMM_WORLD, 0, 1);
  p8est_connectivity_destroy (conn);
  snprintf (file, BUFSIZ, "t8_brick_random_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh, file, 1.);

  t8_cmesh_init (&cmesh_part);
  t8_cmesh_set_partition_from (cmesh_part, cmesh, -1,
                               t8_partition_offset_random (sc_MPI_COMM_WORLD,
                                                           t8_cmesh_get_num_trees
                                                           (cmesh), 0));
  t8_cmesh_commit (cmesh_part);

  if (mpisize > 1 && 1) {
    t8_cmesh_init (&cmesh_part2);

    t8_cmesh_set_partition_from (cmesh_part2, cmesh_part, -1,
                                 t8_partition_offset_random
                                 (sc_MPI_COMM_WORLD,
                                  t8_cmesh_get_num_trees (cmesh), 0));
    t8_cmesh_commit (cmesh_part2);

    snprintf (file, BUFSIZ, "t8_brick_partition_random2_%04d", mpirank);
    t8_cmesh_vtk_write_file (cmesh_part2, file, 1.0);
  }
  else {
    cmesh_part2 = cmesh_part;
  }
  snprintf (file, BUFSIZ, "t8_brick_partition_random_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh_part, file, 1.0);
  t8_cmesh_unref (&cmesh_part2);
}

static void
t8_partition ()
{
  t8_cmesh_t          cmesh, cmesh_part, cmesh_part2;
  char                file[BUFSIZ];
  p4est_connectivity_t *conn;
  int                 mpirank, mpiret, mpisize;

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  conn = p4est_connectivity_new_brick (2, 2, 0, 0);
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, 0, 1);
  p4est_connectivity_destroy (conn);
  snprintf (file, BUFSIZ, "t8_brick_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh, file, 1.);

  t8_cmesh_init (&cmesh_part);
  t8_cmesh_set_partition_from (cmesh_part, cmesh, -1,
                               t8_partition_offset (0, sc_MPI_COMM_WORLD,
                                                    t8_cmesh_get_num_trees
                                                    (cmesh)));
  t8_cmesh_commit (cmesh_part);
  if (mpisize > 1 && 0) {
    t8_cmesh_init (&cmesh_part2);
    t8_cmesh_set_partition_from (cmesh_part2, cmesh_part, -1,
                                 t8_partition_offset (1, sc_MPI_COMM_WORLD,
                                                      t8_cmesh_get_num_trees
                                                      (cmesh)));
    t8_cmesh_commit (cmesh_part2);
    snprintf (file, BUFSIZ, "t8_brick_partition2_%04d", mpirank);
    t8_cmesh_vtk_write_file (cmesh_part2, file, 1.0);
  }
  else {
    cmesh_part2 = cmesh_part;
  }
  snprintf (file, BUFSIZ, "t8_brick_partition_%04d", mpirank);
  t8_cmesh_vtk_write_file (cmesh_part, file, 1.0);
  t8_cmesh_unref (&cmesh_part2);
}

int
main (int argc, char **argv)
{
  int                 mpiret, loop;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  //t8_partition ();
  for (loop = 0; loop < 1; loop++) {
    t8_random_partition ();
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
