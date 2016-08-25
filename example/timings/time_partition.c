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
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <p4est_connectivity.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include "t8_cmesh/t8_cmesh_types.h"

#if 0
/* Repartition a partitioned cmesh by shipping half of the local trees
 * to the next process. */
void t8_time_half (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh_partition;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[1];

  /* Initialize new cmesh and set partitioning */
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_offsets (cmesh_partition,
                                  t8_cmesh_offset_half (cmesh, comm));
  /* Start timer */
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  /* commit (= partition) the new cmesh */
  t8_cmesh_commit (cmesh_partition, comm);
  /* measure passed time */
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "Partition");
  /* print stats */
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);

  /* cleanup */
  t8_cmesh_destroy (&cmesh_partition);
}

/* 1) Take a cmesh from a p4est brick connectivity.
 * 2) Repartition it according to the half partitioning rule.
 * 3) Refine the original cmesh once. level--
 * 4) Goto 2) if level > 0
 *
 * The half partitioning rule means that each process ships half of its
 * trees to the next process. */
void t8_time_brick_refine_half (int x, int y, int x_periodix, int y_periodic,
                                sc_MPI_Comm comm, int level)
{
  t8_cmesh_t    cmesh[2];
  p4est_connectivity_t *p4est_conn;
  int           ref_level, new_ind, old_ind;

  SC_ABORT ("The refine function is not implemented yet.");
  /* Create cmesh from p4est brick connectivity. */
  p4est_conn = p4est_connectivity_new_brick (x, y, x_periodix, y_periodic);
  cmesh[0] = t8_cmesh_new_from_p4est (p4est_conn, comm, 0, 1);
  /* We do not need the p4est connectivity anymore, so we destroy it */
  p4est_connectivity_destroy (p4est_conn);

  /* Time the repartitioning */
  new_ind = 0;
  t8_time_half (cmesh[new_ind], comm);
  /* The refinement and partition loop */
  for (ref_level = 1;ref_level <= level;ref_level++) {
    new_ind = ref_level % 2;
    old_ind = 1 - ref_level % 2;
    /* Initialize the new cmesh that will be the old but refined. */
    t8_cmesh_init (&cmesh[ref_level % 2]);
    /* Set the new cmesh to be partitioned from the old one */
    t8_cmesh_set_derive (cmesh[new_ind], cmesh[old_ind]);
    t8_cmesh_set_refine (cmesh[new_ind], 1);

    /* Unref the old cmesh such that it gets destroyed as soon as the
     * new one is committed. */
    t8_cmesh_unref (&cmesh[old_ind]);
    /* Commit the new, refined cmesh. */
    t8_cmesh_commit (cmesh[new_ind], comm);
    /* vtk output of the refined mesh */
    {
      char filename[BUFSIZ];
      int mpirank, mpiret;

      mpiret = sc_MPI_Comm_rank (comm, &mpirank);
      SC_CHECK_MPI (mpiret);
      snprintf (filename, BUFSIZ, "cmesh_box_partition_level%02i_%04d",
                ref_level, mpirank);
      t8_cmesh_vtk_write_file (cmesh[new_ind], filename, 1.0);
    }
    /* Time the repartitioning */
    t8_time_half (cmesh[new_ind], comm);
  }
  t8_cmesh_destroy (&cmesh[new_ind]);
}
#endif

#if 1
/* Create a cmesh from an x times y p4est box connectivity uniform level 0
 * partitioned. Repartition it shipping half of each processes quadrants to
 * the next process. */
void t8_time_cmesh_partition_brick (int x, int y, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  t8_cmesh_t          cmesh_partition;
  t8_shmem_array_t    new_partition;

  t8_cprofile_t       *profile;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[6];

  /* Create a disjoint brick cmesh with x time y trees on each process */
  cmesh = t8_cmesh_new_disjoint_bricks (x, y, 0, 0, comm);
  /* Allocate profiling struct */
  profile = T8_ALLOC_ZERO (t8_cprofile_t, 1);

  t8_global_productionf ("Committed cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh));
  /* Set up cmesh_partition to be a repartition of cmesh. */
  t8_cmesh_init (&cmesh_partition);
  cmesh_partition->profile = profile;
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  /* The new cmesh is partitioned according to a uniform level 1 refinement */
  new_partition = t8_cmesh_offset_percent (cmesh, comm, 50);
  t8_cmesh_set_partition_offsets (cmesh_partition, new_partition);

  /* Start timer */
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);
  /* commit (= partition) the second cmesh */
  t8_cmesh_commit (cmesh_partition, comm);
  /* measure passed time */
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "Partition");
  t8_global_productionf ("Partitioned cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh_partition));
  sc_stats_set1 (&stats[1], cmesh_partition->profile->partition_trees_shipped,
      "Number of trees sent.");
  sc_stats_set1 (&stats[2], cmesh_partition->profile->partition_ghosts_shipped,
      "Number of ghosts sent.");
  sc_stats_set1 (&stats[3], cmesh_partition->profile->partition_bytes_sent,
      "Number of bytes sent.");
  sc_stats_set1 (&stats[4], cmesh_partition->profile->partition_runtime,
      "Partition runtime (cmesh measured).");
  sc_stats_set1 (&stats[5], cmesh_partition->profile->commit_runtime,
      "Commit runtime (cmesh measured).");
  /* print stats */
  sc_stats_compute (sc_MPI_COMM_WORLD, 6, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 6, stats, 1, 1);
  /* vtk output */
  {
    char filename[BUFSIZ];
    int mpirank, mpiret;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    snprintf (filename, BUFSIZ, "cmesh_box_partition_%04d", mpirank);
    t8_cmesh_vtk_write_file (cmesh_partition, filename, 1.0);
  }
  /* memory clean-up */
  t8_cmesh_destroy (&cmesh);
  t8_cmesh_destroy (&cmesh_partition);
}
#endif

int main (int argc, char *argv[])
{
  int                 mpiret;
  int                 first_argc;
  int                 x_dim, y_dim;
  int                 help = 0;
  int                 refinement_level = 0;
  sc_options_t       *opt;

  /* Initialize MPI, sc, p4est and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Setup for command line options */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'x', "x-dim", &x_dim, 1,
                      "Number of mesh cells in x direction");
  sc_options_add_int (opt, 'y', "y-dim", &y_dim, 1,
                      "Number of mesh cells in y direction");
#if 0
  sc_options_add_int (opt, 'l', "level", &refinement_level, 0,
                      "If > 0 then we refine the coarse mesh as often as l and"
                      " partition on each level.");
#endif
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");

  /* parese command line options */
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                                 opt, argc, argv);
  /* check for wrong usae of arguments */
  if (first_argc < 0 || first_argc != argc
      || x_dim <= 0 || y_dim <= 0 || refinement_level < 0) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    /* Execute this part of the code if all options are correctly set */
    t8_global_productionf ("Starting with x-dim = %i, y-dim = %i\n", x_dim, y_dim);
    t8_time_cmesh_partition_brick (x_dim, y_dim, sc_MPI_COMM_WORLD);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  return 0;
}
