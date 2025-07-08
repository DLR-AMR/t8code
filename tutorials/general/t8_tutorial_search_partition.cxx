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

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_cmesh.h>                           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <sc_options.h>                         /* CLI parser */

typedef struct search_partition_global
{
  /* p4est mesh */
  int uniform_level;  /* level of initial uniform refinement */
  t8_forest_t forest; /* the resulting forest */
} search_partition_global_t;

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  int first_argc, ue;
  sc_options_t *opt;
  search_partition_global_t global, *g = &global;

  /*
   * Init
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  /* Define command line options of this tutorial. */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &g->uniform_level, 3, "Level of uniform refinement");

  /* Proceed in run-once loop for clean abort. */
  ue = 0;
  do {
    /* Parse command line options */
    first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
    if (first_argc < 0) {
      P4EST_GLOBAL_LERROR ("Invalid option format.\n");
      ue = 1;
      break;
    }

    /* Check options for consistency. */
    if (g->uniform_level < 0 || g->uniform_level > 18) {
      /* compared against hard-coded maxlevel of a brick forest */
      P4EST_GLOBAL_LERRORF ("Uniform level out of bounds 0..18\n");
      ue = 1;
    }
    if (ue) {
      break;
    }
    sc_options_print_summary (p4est_get_package_id (), SC_LP_ESSENTIAL, opt);

    /* Print a message on the root process. */
    t8_global_productionf (" [search] \n");
    t8_global_productionf (" [search] Hello, this is the partition search example of t8code.\n");
    t8_global_productionf (
      " [search] We will search for all elements in a forest that contain randomly created particles.\n");
    t8_global_productionf (" [search] \n");

    /*
     *  Build forest and particles.
     */
    comm = sc_MPI_COMM_WORLD;
    /* Build a cube cmesh with tet, hex, and prism trees. */
    cmesh = t8_cmesh_new_brick_3d (2, 2, 2, 0, 0, 0, comm);
    /* Build a uniform forest on it. */
    g->forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), g->uniform_level, 0, comm);

    /* Destroy the forest. */
    t8_forest_unref (&g->forest);
  } while (0);
  if (ue) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  /* Close MPI environment. */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
