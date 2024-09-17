/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

#include <sc_flops.h>
#include <sc_options.h>
#include <sc_statistics.h>

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_eclass.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_helpers.h>

/* In this file benchmarks the `t8_set_join_by_vertices` routine by reading in
 * a given mesh file, retrieving the vertices and building the face
 * connectivity. The benchmark results are compared to the
 * `t8_cmesh_readmshfile` routine. Up to now, the latter is much faster. In the
 * future, `t8_set_join_by_vertices` may be optimized to close up to the reader
 * in terms of speed.
 */

static void
test_with_cmesh (t8_cmesh_t cmesh)
{
  t8_locidx_t ntrees = t8_cmesh_get_num_local_trees (cmesh);

  t8_global_productionf ("ntrees = %d.\n", ntrees);

  /* Arrays for the face connectivity computations via vertices. */
  double *all_verts = T8_ALLOC (double, ntrees *T8_ECLASS_MAX_CORNERS *T8_ECLASS_MAX_DIM);
  t8_eclass_t *all_eclasses = T8_ALLOC (t8_eclass_t, ntrees);

  /* Retrieve all tree vertices and element classes and store them into arrays. */
  for (t8_locidx_t itree = 0; itree < ntrees; itree++) {
    t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    all_eclasses[itree] = eclass;

    const double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);

    const int nverts = t8_eclass_num_vertices[eclass];

    for (int ivert = 0; ivert < nverts; ivert++) {
      for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
        all_verts[T8_3D_TO_1D (ntrees, T8_ECLASS_MAX_CORNERS, T8_ECLASS_MAX_DIM, itree, ivert, icoord)]
          = vertices[T8_2D_TO_1D (nverts, T8_ECLASS_MAX_DIM, ivert, icoord)];
      }
    }
  }

  sc_flopinfo_t fi, snapshot;
  sc_statinfo_t stats[1];

  /* Start timer */
  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  /* Compute face connectivity. */
  const int do_both_directions = 0;
  t8_cmesh_set_join_by_vertices (NULL, ntrees, all_eclasses, all_verts, NULL, do_both_directions);

  /* Measure passed time. */
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "t8_cmesh_set_join_by_vertices");

  /* Print stats. */
  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);

  T8_FREE (all_verts);
  T8_FREE (all_eclasses);
}

int
main (int argc, char **argv)
{
  char usage[BUFSIZ];
  /* brief help message */
  int sreturnA = snprintf (usage, BUFSIZ,
                           "Usage:\t%s <OPTIONS>\n\t%s -h\t"
                           "for a brief overview of all options.",
                           basename (argv[0]), basename (argv[0]));

  char help[BUFSIZ];
  /* long help message */
  int sreturnB = snprintf (help, BUFSIZ, "Profile `t8_cmesh_set_join_by_vertices` via given mesh file.\n\n%s\n", usage);

  if (sreturnA > BUFSIZ || sreturnB > BUFSIZ) {
    /* The usage string or help message was truncated */
    /* Note: gcc >= 7.1 prints a warning if we 
     * do not check the return value of snprintf. */
    t8_debugf ("Warning: Truncated usage string and help message to '%s' and '%s'\n", usage, help);
  }

  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  int helpme;

  const char *meshfile;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_string (opt, 'f', "fileprefix", &meshfile, NULL, "File prefix of the mesh file (without .msh)");

  int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && meshfile != NULL) {
    t8_global_productionf ("meshfile = %s\n", meshfile);

    const int dim = 3;
    const int main_proc = 0;
    const int partition = 0;
    const int use_cad_geometry = 0;

    t8_cmesh_t cmesh;

    {
      sc_flopinfo_t fi, snapshot;
      sc_statinfo_t stats[1];

      /* Start timer */
      sc_flops_start (&fi);
      sc_flops_snap (&fi, &snapshot);

      cmesh = t8_cmesh_from_msh_file (meshfile, partition, sc_MPI_COMM_WORLD, dim, main_proc, use_cad_geometry);

      /* Measure passed time. */
      sc_flops_shot (&fi, &snapshot);
      sc_stats_set1 (&stats[0], snapshot.iwtime, "t8_cmesh_from_msh_file");
      /* Print stats. */
      sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
      sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);
    }

    test_with_cmesh (cmesh);

    t8_cmesh_unref (&cmesh);
  }
  else {
    /* Display help message and usage. */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
