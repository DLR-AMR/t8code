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
#include <sc_functions.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <p4est_connectivity.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_vtk/t8_vtk_writer.h>

#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh/t8_cmesh_cad.hxx>
#include <t8_cmesh_readmshfile.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <example/common/t8_example_common.hxx>
#include <t8_types/t8_vec.hxx>

/* This is the user defined data used to define the
 * region in which we partition.
 * Let E be a plane through the origin, then we define all cells
 * in between the two planes c_min*E and c_max*E.
 */
typedef struct
{
  double c_min, c_max; /* constants that define the thickness of the refinement region */
  t8_3D_vec normal;    /* normal vector to the plane E */
  int base_level;      /* A given level that is not coarsend further, see -l argument */
  int max_level;       /* A max level that is not refined further, see -L argument */
} adapt_data_t;

/* refine the forest in a band, given by a plane E and two constants
 * c_min, c_max. We refine the cells in the band c_min*E, c_max*E */
static int
t8_band_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_eclass_t tree_class,
               [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  int level, base_level, max_level;
  t8_3D_vec elem_midpoint;
  adapt_data_t *adapt_data;

  T8_ASSERT (!is_family || num_elements == scheme->element_get_num_children (tree_class, elements[0]));
  level = scheme->element_get_level (tree_class, elements[0]);
  /* Get the minimum and maximum x-coordinate from the user data pointer of forest */
  adapt_data = (adapt_data_t *) t8_forest_get_user_data (forest);
  t8_3D_vec normal = adapt_data->normal;
  base_level = adapt_data->base_level;
  max_level = adapt_data->max_level;
  /* Compute the coordinates of the anchor node. */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint.data ());

  /* Calculate elem_midpoint - c_min n */
  t8_axpy (normal, elem_midpoint, -adapt_data->c_min);

  /* The purpose of the factor C*h is that the levels get smaller, the
   * closer we get to the interface. We refine a cell if it is at most
   * C times its own height away from the interface */
  if (t8_dot (elem_midpoint, normal) >= 0) {
    /* if the anchor node is to the right of c_min*E,
     * check if it is to the left of c_max*E */

    /* set elem_midpoint to the original anchor - c_max*normal */
    t8_axpy (normal, elem_midpoint, adapt_data->c_min - adapt_data->c_max);
    if (t8_dot (elem_midpoint, normal) <= 0) {
      if (level < max_level) {
        /* We do refine if level smaller 1+base level and the anchor is
         * to the left of c_max*E */
        return 1;
      }
    }
    else if (is_family && level > base_level) {
      /* Otherwise, we coarse if we have a family and level is greater
       * than the base level. */
      return -1;
    }
  }
  else if (is_family && level > base_level) {
    /* If element lies out of the refinement region and a family was given
     * as argument, we coarsen to level base level */
    /* set elem_midpoint to the original midpoint - c_max*normal */
    return -1;
  }
  return 0;
}

/* Create a cmesh from a .msh files uniform level 0
 * partitioned. */
static void
t8_time_forest_cmesh_mshfile (t8_cmesh_t cmesh, const char *vtu_prefix, sc_MPI_Comm comm, int init_level, int max_level,
                              int no_vtk, double x_min_max[2], double T, double delta_t, int do_ghost, int do_balance)
{
  char forest_vtu[BUFSIZ], cmesh_vtu[BUFSIZ];
  adapt_data_t adapt_data;
  t8_forest_t forest, forest_adapt, forest_partition;
  double t;
  int cmesh_is_partitioned;
  int time_step;
  double adapt_time = 0, ghost_time = 0, partition_time = 0, new_time = 0, total_time = 0, balance_time = 0;
  sc_statinfo_t times[6];

  t8_global_productionf ("Committed cmesh with %lli global trees.\n", (long long) t8_cmesh_get_num_trees (cmesh));

  /* If the input cmesh is partitioned then we use a partitioned cmehs
   * and also repartition it in each timestep (happens automatically in
   * t8_forest_commit). We have to initially start with a uniformly refined
   * cmesh in order to be able to construct the forest on it.
   * If on the other hand, the input cmesh was replicated, then we keep it
   * as replicated throughout. */
  sc_stats_init (&times[0], "new");
  sc_stats_init (&times[1], "adapt");
  sc_stats_init (&times[2], "ghost");
  sc_stats_init (&times[3], "partition");
  sc_stats_init (&times[4], "balance");
  sc_stats_init (&times[5], "total");
  total_time -= sc_MPI_Wtime ();

  cmesh_is_partitioned = t8_cmesh_is_partitioned (cmesh);
  /* Initialize forest and set cmesh */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh, comm);
  /* Set the element scheme */
  t8_forest_set_scheme (forest, t8_scheme_new_default ());
  /* Set the initial refinement level */
  t8_forest_set_level (forest, init_level);
  /* Commit the forest */
  new_time -= sc_MPI_Wtime ();
  t8_forest_commit (forest);
  new_time += sc_MPI_Wtime ();
  sc_stats_set1 (&times[0], new_time, "new");
  /* Set the permanent data for adapt. */
  adapt_data.normal[0] = 0.8;
  adapt_data.normal[1] = 0.3;
  adapt_data.normal[2] = 0.0;
  t8_normalize (adapt_data.normal);
  adapt_data.base_level = init_level;
  adapt_data.max_level = max_level;
  /* Start the time loop, in each time step the refinement front  moves
   * further through the domain */
  for (t = 0, time_step = 0; t < T; t += delta_t, time_step++) {
    /* Adapt the forest */
    /* TODO: profiling */
    t8_forest_init (&forest_adapt);
    t8_forest_set_adapt (forest_adapt, forest, t8_band_adapt, 1);
    t8_forest_set_profiling (forest_adapt, 1);
    /* Set the minimum and maximum x-coordinates as user data */
    adapt_data.c_min = x_min_max[0] + t;
    adapt_data.c_max = x_min_max[1] + t;
    t8_forest_set_user_data (forest_adapt, (void *) &adapt_data);
    t8_forest_commit (forest_adapt);
    t8_forest_compute_profile (forest_adapt);
    t8_forest_ref (forest_adapt);

    /* partition the adapted forest */
    t8_forest_init (&forest_partition);
    /* partition the adapted forest */
    t8_forest_set_partition (forest_partition, forest_adapt, 0);

    /* If desired, create ghost elements and balance */
    t8_forest_set_profiling (forest_partition, 1);
    if (do_ghost) {
      t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
    }
    if (do_balance) {
      t8_forest_set_balance (forest_partition, NULL, 0);
    }
    t8_forest_commit (forest_partition);
    t8_forest_compute_profile (forest_partition);
    t8_cmesh_print_profile (t8_forest_get_cmesh (forest_partition));
    /* Set forest to the partitioned forest, so it gets adapted
     * in the next time step. */
    forest = forest_partition;

    /* Set the vtu output name */
    if (!no_vtk) {
      snprintf (forest_vtu, BUFSIZ, "%s_forest_partition_%03d", vtu_prefix, time_step);
      snprintf (cmesh_vtu, BUFSIZ, "%s_cmesh_partition_%03d", vtu_prefix, time_step);
      t8_forest_write_vtk (forest_partition, forest_vtu);
      t8_cmesh_vtk_write_file (t8_forest_get_cmesh (forest_partition), cmesh_vtu);
      t8_debugf ("Wrote partitioned forest and cmesh\n");
    }
    if (cmesh_is_partitioned) {
      /* Print runtimes and statistics of forest and cmesh partition */
      t8_cmesh_print_profile (t8_forest_get_cmesh (forest_partition));
    }
    t8_forest_print_profile (forest_partition);
    /* clean-up */
    t8_forest_unref (&forest_adapt);
    /* TIME-LOOP ends here */
  }
  /* memory clean-up */
  total_time += sc_MPI_Wtime ();
  sc_stats_accumulate (&times[0], new_time);
  sc_stats_accumulate (&times[1], adapt_time);
  sc_stats_accumulate (&times[2], ghost_time);
  sc_stats_accumulate (&times[3], partition_time);
  sc_stats_accumulate (&times[4], balance_time);
  sc_stats_accumulate (&times[5], total_time);
  sc_stats_compute (comm, 6, times);
  sc_stats_print (t8_get_package_id (), SC_LP_ESSENTIAL, 6, times, 1, 1);
  t8_forest_unref (&forest_partition);
}

#undef USE_CMESH_PARTITION

/* Construct a cmesh either from a .msh mesh file or from a
 * collection of cmesh files constructed with t8_cmesh_save.
 * If msh_file is NULL, the cmesh is loaded from the cmesh_file and num_files
 * must be specified. If cmesh_file is NULL, the cmesh is loaded from the .msh
 * file and mesh_dim must be specified. */
t8_cmesh_t
t8_time_forest_create_cmesh (const char *msh_file, int mesh_dim, const char *cmesh_file, int num_files,
                             sc_MPI_Comm comm, int init_level, int stride, int use_cad)
{
  t8_cmesh_t cmesh;
  t8_cmesh_t cmesh_partition;
  int partition;

  T8_ASSERT (msh_file == NULL || cmesh_file == NULL);

  if (msh_file != NULL) {
    if (use_cad) {
      partition = 0;
      t8_global_productionf ("The cmesh is not partitioned due to the usage of the curved mesh option. \n"
                             "Timing will not be comparable to non-curved meshes. \n");
    }
    else {
      partition = 1;
    }
    /* Create a cmesh from the given mesh files */
    cmesh = t8_cmesh_from_msh_file ((char *) msh_file, partition, comm, mesh_dim, 0, use_cad);
  }
  else {
    T8_ASSERT (cmesh_file != NULL);
    SC_CHECK_ABORT (num_files > 0, "Must specify valid number of files.\n");
    /* Load the cmesh from the stored files and evenly distribute it
     * among all ranks */
    cmesh = t8_cmesh_load_and_distribute (cmesh_file, num_files, comm, T8_LOAD_STRIDE, stride);
    /* Partition only if more than 1 input file */
    partition = num_files > 1;
  }
  SC_CHECK_ABORT (cmesh != NULL, "Error when creating cmesh.\n");

  if (partition) {
    /* partition the cmesh uniformly */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_partition, init_level, t8_scheme_new_default ());
    t8_cmesh_set_profiling (cmesh_partition, 1);
    t8_cmesh_commit (cmesh_partition, comm);
    return cmesh_partition;
  }
  return cmesh;
}

int
main (int argc, char *argv[])
{
  int mpiret, mpisize;
  int first_argc;
  int level, level_diff;
  int help = 0, no_vtk, do_ghost, do_balance, use_cad;
  int dim, num_files;
  int test_tet, test_linear_cylinder, test_cad_cylinder, test_hybrid_cube, test_hex_cube;
  int stride;
  int cmesh_level;
  double T, delta_t, cfl;
  sc_options_t *opt;
  t8_cmesh_t cmesh;
  const char *mshfileprefix, *cmeshfileprefix;
  const char *vtu_prefix;
  double x_min_max[2]; /* At position 0 the minimum x-coordinate,
                                           at position 1 the maximum x-coordinate. */

  /* Initialize MPI, sc, p4est and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* get mpisize */
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_STATISTICS);

  /* Setup for command line options */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &help, "Display a short help message.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk, "Do not write vtk output.");
  sc_options_add_string (opt, 'f', "mshfile", &mshfileprefix, NULL,
                         "If specified, the cmesh is constructed from a .msh file with the given prefix. "
                         "The files must end in .msh and be created with gmsh.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2, "Together with -f: The dimension of the coarse mesh. 2 or 3.");
  sc_options_add_string (opt, 'c', "cmeshfile", &cmeshfileprefix, NULL,
                         "If specified, the cmesh is constructed from a collection of cmesh files. "
                         "Created with t8_cmesh_save. The number of files must then be specified with the -n option.");
  sc_options_add_int (opt, 'n', "nfiles", &num_files, -1,
                      "If the -c option is used, the number of cmesh files must be specified as an argument here. "
                      "If n=1 then the cmesh will be replicated throughout the test.");
  sc_options_add_int (opt, 's', "stride", &stride, 16,
                      "If -c and -n are used, only every s-th MPI rank will read a .cmesh file (file number: rank/s)."
                      "Default is 16.");
  sc_options_add_switch (opt, 't', "test-tet", &test_tet,
                         "Use a cmesh that tests all tet face-to-face connections."
                         " If this option is used -o is enabled automatically. Not allowed with -f and -c.");
  sc_options_add_switch (opt, 'L', "test-linear-cylinder", &test_linear_cylinder,
                         "Use a linear cmesh to compare linear and cad geometry performance."
                         " If this option is used -o is enabled automatically. Not allowed with -f and -c.");
  sc_options_add_switch (opt, 'O', "test-cad-cylinder", &test_cad_cylinder,
                         "Use a cad cmesh to compare linear and cad geometry performance."
                         " If this option is used -o is enabled automatically. Not allowed with -f and -c.");
  sc_options_add_switch (opt, 'C', "test-hybridcube", &test_hybrid_cube,
                         "Use a hypercube with Tet, Prism and Hex elements as cmesh."
                         " If this option is used -o is enabled automatically. Not allowed with -f and -c.");
  sc_options_add_switch (opt, 'H', "test-hexcube", &test_hex_cube,
                         "Use a hypercube with Hex elements as cmesh."
                         " If this option is used -o is enabled automatically. Not allowed with -f and -c.");
  sc_options_add_int (opt, 'l', "level", &level, 0, "The initial uniform refinement level of the forest.");
  sc_options_add_int (opt, 'r', "rlevel", &level_diff, 1,
                      "The number of levels that the forest is refined from the initial level.");
  sc_options_add_int (opt, '\0', "cmesh-level", &cmesh_level, -1,
                      "Starting level of the linear or cad cmesh, default is 0. Only viable with -L or -O.");
  sc_options_add_double (opt, 'x', "xmin", x_min_max, 0, "The minimum x coordinate in the mesh.");
  sc_options_add_double (opt, 'X', "xmax", x_min_max + 1, 1, "The maximum x coordinate in the mesh.");
  sc_options_add_double (opt, 'T', "time", &T, 1,
                         "The simulated time span. We simulate the time from 0 to T. T has to be > 0.");
  sc_options_add_double (opt, 'D', "delta_t", &delta_t, 0.08,
                         "The time step in each simulation step. Deprecated, use -C instead.");
  /* CFL number. delta_t = CFL * 0.64 / 2^level */
  sc_options_add_double (opt, 'C', "cfl", &cfl, 0,
                         "The CFL number. If specified, then delta_t is set to CFL * 0.64 / 2^level. "
                         "Overwrites any other delta_t setting.");
  sc_options_add_switch (opt, 'g', "ghost", &do_ghost, "Create ghost elements.");
  sc_options_add_switch (opt, 'b', "balance", &do_balance, "Establish a 2:1 balance in the forest.");
  sc_options_add_switch (opt, 'z', "use_cad", &use_cad,
                         "If used, meshes will be curved to original geometries (msh- and brep-files necessary).");

  /* parse command line options */
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
  /* check for wrong usage of arguments */
  if (first_argc < 0 || first_argc != argc || dim < 2 || dim > 3
      || (cmeshfileprefix == NULL && mshfileprefix == NULL && test_tet == 0 && test_cad_cylinder == 0
          && test_linear_cylinder == 0 && test_hybrid_cube == 0 && test_hex_cube == 0)
      || stride <= 0 || (num_files - 1) * stride >= mpisize || cfl < 0 || T <= 0
      || test_tet + test_linear_cylinder + test_cad_cylinder + test_hybrid_cube + test_hex_cube > 1
      || (cmesh_level >= 0 && (!test_linear_cylinder && !test_cad_cylinder && !test_hybrid_cube && !test_hex_cube))
      || ((mshfileprefix != NULL || cmeshfileprefix != NULL)
          && (test_linear_cylinder || test_cad_cylinder || test_tet || test_hybrid_cube || test_hex_cube))
      || (mshfileprefix == NULL && use_cad)) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    /* Execute this part of the code if all options are correctly set */
    /* Set correct timestep, if -C was specified with a value greater 0, then
     * we overwrite the delta_t setting. We choose this to be backwards compatible to
     * call that use the delta_t option. Eventually, we will remove the delta_t option
     * completely. */
    if (cfl > 0) {
      delta_t = cfl * 0.64 / (1 << level);
    }
    t8_global_productionf ("Using delta_t = %f\n", delta_t);
    if (mshfileprefix != NULL) {
      cmesh = t8_time_forest_create_cmesh (mshfileprefix, dim, NULL, -1, sc_MPI_COMM_WORLD, level, stride, use_cad);
      vtu_prefix = mshfileprefix;
    }
    else if (test_tet) {
      cmesh = t8_cmesh_new_tet_orientation_test (sc_MPI_COMM_WORLD);
      vtu_prefix = "test_tet";
    }
    else if (test_linear_cylinder || test_cad_cylinder) {
      if (cmesh_level < 0) {
        cmesh_level = 0;
      }
      cmesh = t8_cmesh_new_hollow_cylinder (sc_MPI_COMM_WORLD, 4 * sc_intpow (2, cmesh_level),
                                            sc_intpow (2, cmesh_level), sc_intpow (2, cmesh_level), test_cad_cylinder);
      test_linear_cylinder ? vtu_prefix = "test_linear_cylinder" : vtu_prefix = "test_cad_cylinder";
    }
    else if (test_hybrid_cube) {
      cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
      vtu_prefix = "test_hypercube_hybrid";
    }
    else if (test_hex_cube) {
      cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
      vtu_prefix = "test_hypercube_hex";
    }
    else {
      T8_ASSERT (cmeshfileprefix != NULL);
      cmesh
        = t8_time_forest_create_cmesh (NULL, -1, cmeshfileprefix, num_files, sc_MPI_COMM_WORLD, level, stride, use_cad);
      vtu_prefix = cmeshfileprefix;
    }
    t8_time_forest_cmesh_mshfile (cmesh, vtu_prefix, sc_MPI_COMM_WORLD, level, level + level_diff, no_vtk, x_min_max, T,
                                  delta_t, do_ghost, do_balance);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
