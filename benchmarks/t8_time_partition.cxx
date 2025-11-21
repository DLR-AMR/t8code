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
#include <t8_cmesh/t8_cmesh.h>
#include <t8_vtk/t8_vtk_writer.h>

#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* Repartition a partitioned cmesh by shipping half of the local trees
 * to the next process. */
void
t8_time_half (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh_partition;
  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  /* Initialize new cmesh and set partitioning */
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_set_partition_offsets (cmesh_partition, t8_cmesh_offset_half (cmesh, comm));
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

#if 1
/* add x*mpirank to all of the x-coordinates of a tree in cmesh.
 * This is useful to separate the disjoint parts visually. */
void
t8_time_cmesh_translate_coordinates (t8_cmesh_t cmesh, double x, sc_MPI_Comm comm)
{
  double *vertices;
  t8_locidx_t itree;
  int ivertex, mpiret, mpirank;
  t8_eclass_t eclass;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  for (itree = 0; itree < t8_cmesh_get_num_local_trees (cmesh); itree++) {
    /* loop over all trees and get a pointer to their vertices */
    vertices = (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), T8_CMESH_VERTICES_ATTRIBUTE_KEY, itree);
    eclass = t8_cmesh_get_tree_class (cmesh, itree);
    for (ivertex = 0; ivertex < t8_eclass_num_vertices[eclass]; ivertex++) {
      /* For each tree vertex, translate its x-coordinate */
      vertices[3 * ivertex] += x * mpirank;
    }
  }
}

/* Create a cmesh from an x times y p4est box connectivity uniform level 0
 * partitioned. Repartition it by shipping 43% of each processes quadrants to
 * the next process. */
void
t8_time_cmesh_partition_brick (int x, int y, int z, sc_MPI_Comm comm, int no_vtk)
{
  t8_cmesh_t cmesh;
  t8_cmesh_t cmesh_partition;
  t8_shmem_array_t new_partition;

  t8_cprofile_t *profile;

  /* Create a disjoint brick cmesh with x time y trees on each process */
  cmesh = t8_cmesh_new_disjoint_bricks (x, y, z, 1, 1, 1, comm);
  /* Allocate profiling struct */
  profile = T8_ALLOC_ZERO (t8_cprofile_t, 1);

  t8_global_productionf ("Committed cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh));
  /* If output is active, translate the disjoint parts to make them
   * physically disjoint and distinguishable */
  if (!no_vtk) {
    t8_time_cmesh_translate_coordinates (cmesh, x + x / 4., comm);
  }
  /* Set up cmesh_partition to be a repartition of cmesh. */
  t8_cmesh_init (&cmesh_partition);
  cmesh_partition->profile = profile;
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  /* The new cmesh is partitioned according to a uniform level 1 refinement */
  new_partition = t8_cmesh_offset_percent (cmesh, comm, 43);
  t8_cmesh_set_partition_offsets (cmesh_partition, new_partition);
  /* activate profilling for cmesh to obtain run times */
  t8_cmesh_set_profiling (cmesh_partition, 1);
  /* commit (= partition) the second cmesh */
  t8_cmesh_commit (cmesh_partition, comm);
  t8_global_productionf ("Partitioned cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh_partition));
  /* Print run times and statistics */
  t8_cmesh_print_profile (cmesh_partition);

  /* vtk output */
  if (!no_vtk) {
    int mpirank, mpiret;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    t8_cmesh_vtk_write_file (cmesh_partition, "cmesh_box_partition");
  }
  /* memory clean-up */
  t8_cmesh_destroy (&cmesh_partition);
}
#endif

int
main (int argc, char *argv[])
{
  int mpiret;
  int first_argc;
  int x_dim, y_dim, z_dim;
  int no_vtk = 0;
  int help = 0;
  int dim;
  sc_options_t *opt;

  /* Initialize MPI, sc, p4est and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* Setup for command line options */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'x', "x-dim", &x_dim, 1, "Number of mesh cells in x direction.");
  sc_options_add_int (opt, 'y', "y-dim", &y_dim, 1, "Number of mesh cells in y direction.");
  sc_options_add_int (opt, 'z', "z-dim", &z_dim, 0,
                      "Number of mesh cells in z direction."
                      " If specified, then the mesh is automatically 3d.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2,
                      "The dimension of the coarse mesh."
                      " 2 for a quadmesh (z is ignored) and 3 for a hexmesh.");
  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "Disable vtk output");

  /* parse command line options */
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  /* check for wrong usage of arguments */
  if (first_argc < 0 || first_argc != argc || x_dim <= 0 || y_dim <= 0 || dim < 2 || dim > 3) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    if (z_dim <= 0) {
      /* For 2d problems, there are no trees in z direction. */
      dim = 2;
      z_dim = 0;
    }
    else {
      dim = 3;
    }
    /* Execute this part of the code if all options are correctly set */
    t8_global_productionf ("Starting with x-dim = %i, y-dim = %i z-dim = %i\n", x_dim, y_dim, z_dim);
    t8_time_cmesh_partition_brick (x_dim, y_dim, z_dim, sc_MPI_COMM_WORLD, no_vtk);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  return 0;
}
