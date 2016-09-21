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
#include <t8_cmesh_readmshfile.h>
#include <t8_forest.h>
#include <t8_default.h>

double              x_min, x_max;

static void
t8_anchor_element (t8_forest_t forest,
                   t8_eclass_scheme_t * ts, t8_element_t * element,
                   double elem_anchor_f[3])
{
  int                 elem_anchor[3], maxlevel, i;

  /* get the element anchor node */
  t8_element_anchor (ts, element, elem_anchor);
  maxlevel = t8_element_maxlevel (ts);
  for (i = 0; i < 3; i++) {
    /* Calculate the anchor coordinate in [0,1]^3 */
    elem_anchor_f[i] = elem_anchor[i] / (1 << maxlevel);
  }
}

static int
t8_basic_adapt (t8_forest_t forest, t8_topidx_t which_tree,
                t8_eclass_scheme_t * ts,
                int num_elements, t8_element_t * elements[])
{
  int                 level;
  double              elem_anchor[3];

  T8_ASSERT (num_elements == 1 || num_elements ==
             t8_eclass_num_children[ts->eclass]);
  level = t8_element_level (ts, elements[0]);

  t8_anchor_element (forest, ts, elements[0], elem_anchor);

  if (level < 3 && elem_anchor[0] >= x_min && elem_anchor[0] <= x_max) {
    /* refine if level is smaller 3 and
     * anchor x coordinate is in correct range */
    return 1;
  }
  return 0;
}

#if 1
#define NUM_STATS 9
#define USE_CMESH_PARTITION 1   /* Set this to false to use a replicated cmesh
                                 * cmesh throughout the function */
/* Create a cmesh from a .msh files uniform level 0
 * partitioned. */
void
t8_time_forest_cmesh_mshfile (const char *msh_file, int mesh_dim,
                              sc_MPI_Comm comm, int init_level, int no_vtk)
{
  t8_cmesh_t          cmesh;
  t8_cmesh_t          cmesh_partition;
  char                forest_vtu[BUFSIZ], cmesh_vtu[BUFSIZ];

  t8_forest_t         forest, forest_adapt, forest_partition;

  /* Create a cmesh from the given mesh files */
  cmesh = t8_cmesh_from_msh_file ((char *) msh_file, 1, comm, mesh_dim, 0);

  t8_global_productionf ("Committed cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh));
#if USE_CMESH_PARTITION
  /* Set up cmesh_partition to be a repartition of cmesh. */
  t8_cmesh_init (&cmesh_partition);
  t8_cmesh_set_derive (cmesh_partition, cmesh);
  t8_cmesh_unref (&cmesh);
  /* The new cmesh is partitioned according to a uniform init_level refinement */
  t8_cmesh_set_partition_uniform (cmesh_partition, init_level);
  t8_cmesh_commit (cmesh_partition, comm);
#else
  cmesh_partition = cmesh;
#endif
  /* Initialize forest and set cmesh */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh_partition, comm);
  /* Set the element scheme */
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  /* Set the initial refinement level */
  t8_forest_set_level (forest, init_level);
  /* Commit the forest */
  t8_forest_commit (forest);
  /* Adapt the forest */
  t8_forest_init (&forest_adapt);
  t8_forest_set_adapt (forest_adapt, forest, t8_basic_adapt, NULL, 1);
  /* Commit the adapted forest */
  t8_forest_commit (forest_adapt);
  /* write vtk files for adapted forest and cmesh */
  if (!no_vtk) {
    snprintf (forest_vtu, BUFSIZ, "%s_forest_adapt", msh_file);
    snprintf (cmesh_vtu, BUFSIZ, "%s_cmesh_adapt", msh_file);
    t8_forest_write_vtk (forest_adapt, forest_vtu);
    t8_cmesh_vtk_write_file (t8_forest_get_cmesh (forest_adapt),
                             cmesh_vtu, 1.0);
  }
  /* partition the adapted forest */
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest_adapt, 0);
  /* enable profiling for the partitioned forest */
  t8_forest_set_profiling (forest_partition, 1);
  t8_forest_commit (forest_partition);
#if USE_CMESH_PARTITION
  /* Repartition the cmesh of the forest */
  t8_forest_partition_cmesh (forest_partition, comm, 1);
#endif
  /* Set the vtu output name */
  if (!no_vtk) {
    snprintf (forest_vtu, BUFSIZ, "%s_forest_partition", msh_file);
    snprintf (cmesh_vtu, BUFSIZ, "%s_cmesh_partition", msh_file);
    t8_forest_write_vtk (forest_partition, forest_vtu);
    t8_cmesh_vtk_write_file (t8_forest_get_cmesh (forest_partition),
                             cmesh_vtu, 1.0);
  }
  /* Set the vtu output name */
  snprintf (forest_vtu, BUFSIZ, "%s_forest_partition", msh_file);
  snprintf (cmesh_vtu, BUFSIZ, "%s_cmesh_partition", msh_file);
  t8_debugf ("Wrote partitioned forest and cmesh\n");
  t8_forest_write_vtk (forest_partition, forest_vtu);
  t8_cmesh_vtk_write_file (t8_forest_get_cmesh (forest_partition), cmesh_vtu,
                           1.0);

  /* Print runtimes and statistics of forest and cmesh partition */
  t8_forest_print_profile (forest_partition);
  t8_cmesh_print_profile (t8_forest_get_cmesh (forest_partition));
  /* memory clean-up */
  t8_forest_unref (&forest_partition);
  t8_cmesh_destroy (&cmesh_partition);
}

#undef NUM_STATS
#undef USE_CMESH_PARTITION
#endif

int
main (int argc, char *argv[])
{
  int                 mpiret;
  int                 first_argc;
  int                 level;
  int                 help = 0, no_vtk;
  int                 dim;
  sc_options_t       *opt;
  const char         *fileprefix;

  /* Initialize MPI, sc, p4est and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_INFO);

  /* Setup for command line options */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'd', "dim", &dim, 2,
                      "The dimension of the coarse mesh. 2 or 3.");
  sc_options_add_switch (opt, 'h', "help", &help,
                         "Display a short help message.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Do not write vtk output.");
  sc_options_add_string (opt, 'f', "file", &fileprefix, NULL,
                         "The input mesh "
                         "file prefix. The files must end in .msh and be "
                         "created with gmsh.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The initial uniform "
                      "refinement level of the forest.");
  sc_options_add_double (opt, 'x', "xmin", &x_min, 0,
                         "The minimum x coordinate " "in the mesh.");
  sc_options_add_double (opt, 'X', "xmax", &x_max, 1,
                         "The maximum x coordinate " "in the mesh.");

  /* parse command line options */
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                                 opt, argc, argv);
  /* check for wrong usage of arguments */
  if (first_argc < 0 || first_argc != argc || dim < 2 || dim > 3) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    /* Execute this part of the code if all options are correctly set */
    t8_time_forest_cmesh_mshfile (fileprefix, dim, sc_MPI_COMM_WORLD, level,
                                  no_vtk);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  return 0;
}
