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
#include <t8_forest/t8_forest_general.h>        /* Forest definition and basic interface. */
#include <t8_forest/t8_forest_geometrical.h>    /* Element center computation. */
#include <t8_cmesh.h>                           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <sc_options.h>                         /* CLI parser */

typedef struct search_partition_global
{
  /* forest mesh */
  double a[3], b[3], c[3]; /* refinement centers */
  int uniform_level;       /* level of initial uniform refinement */
  int max_level;           /* maximum level of adaptive refinement */
  t8_forest_t forest;      /* the resulting forest */

  /* MPI */
  int mpirank; /* the processes rank */
} search_partition_global_t;

/* Map coordinates from the reference coordinate system of the 2x2x2 brick
 * forest to the unit cube. */
static void
map_coordinates (double *xyz, t8_locidx_t which_tree, double *mapped_xyz)
{
  T8_ASSERT (which_tree < 8); /* assert we have a 2x2(x2) brick */
  mapped_xyz[0] = 0.5 * xyz[0];
  mapped_xyz[1] = 0.5 * xyz[1];
  mapped_xyz[2] = 0.5 * xyz[2];
}

int
t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                   [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                   [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                   [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  double center[3], mapped_center[3];
  double dist, min_dist;

  search_partition_global_t *g = (search_partition_global_t *) t8_forest_get_user_data (forest);

  /* Compute the element center's position in the unit cube. */
  T8_ASSERT (num_elements == 1);
  t8_forest_element_centroid (forest, which_tree, elements[0], center);
  map_coordinates (center, which_tree, mapped_center);

  /* Compute distance to point a. */
  dist = (g->a[0] - mapped_center[0]) * (g->a[0] - mapped_center[0])
         + (g->a[1] - mapped_center[1]) * (g->a[1] - mapped_center[1])
         + (g->a[2] - mapped_center[2]) * (g->a[2] - mapped_center[2]);
  min_dist = sqrt (dist);

  /* Compute distance to point b. */
  dist = (g->b[0] - mapped_center[0]) * (g->b[0] - mapped_center[0])
         + (g->b[1] - mapped_center[1]) * (g->b[1] - mapped_center[1])
         + (g->b[2] - mapped_center[2]) * (g->b[2] - mapped_center[2]);
  min_dist = SC_MIN (min_dist, sqrt (dist));

  /* refine if quadrant center is close enough to either point a or point b */
  return (scheme->element_get_level (tree_class, elements[0])
          < g->max_level - floor (min_dist * (g->max_level - g->uniform_level) / 0.2));
}

static void
create_forest (search_partition_global_t *g)
{
  int mpiret;
  int il;
  t8_gloidx_t old_gnl;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh;
  t8_forest_t forest_adapt;

  /* Get the MPI rank. */
  comm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (comm, &g->mpirank);
  SC_CHECK_MPI (mpiret);

  /* Build a 2x2x2 cube cmesh. */
  cmesh = t8_cmesh_new_brick_3d (2, 2, 2, 0, 0, 0, comm);
  /* Build a uniform forest on it. */
  g->forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), g->uniform_level, 0, comm);

  /* Refine the forest around two refinement centers. */
  for (il = g->uniform_level; il < g->max_level; il++) {
    /* Store global num leaves for future comparison. */
    old_gnl = t8_forest_get_global_num_leaf_elements (g->forest);

    /* Adapt the forest. */
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, g);
    t8_forest_set_adapt (forest_adapt, g->forest, t8_adapt_callback, 0);
    t8_forest_set_partition (forest_adapt, g->forest, 0);
    t8_forest_commit (forest_adapt);
    g->forest = forest_adapt;

    /* Leave loop if no refinement occured. */
    if (old_gnl == t8_forest_get_global_num_leaf_elements (g->forest)) {
      break;
    }
  }
}

static void
cleanup (search_partition_global_t *g)
{
  /* Destroy the forest. */
  t8_forest_unref (&g->forest);
}

static void
run (search_partition_global_t *g)
{
  /* Create a 2x2x2 brick forest covering the unit square. */
  create_forest (g);

  /* Fee memory. */
  cleanup (g);
}

int
main (int argc, char **argv)
{
  int mpiret;
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
  sc_options_add_int (opt, 'L', "maxlevel", &g->max_level, 5, "Level of maximum refinement");

  /* Proceed in run-once loop for clean abort. */
  ue = 0;
  do {
    /* Parse command line options */
    first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);
    if (first_argc < 0) {
      t8_global_errorf ("Invalid option format.\n");
      ue = 1;
      break;
    }

    /* Check options for consistency. */
    if (g->uniform_level < 0 || g->uniform_level > 18) {
      /* compared against hard-coded maxlevel of a brick forest */
      t8_global_errorf ("Uniform level out of bounds 0..18\n");
      ue = 1;
    }
    if (g->max_level < 0 || g->max_level > 18) {
      t8_global_errorf ("Maximum level out of bounds 0..18\n");
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
      " [search] We will search for all elements in a forest that contains randomly created particles.\n");
    t8_global_productionf (" [search] \n");

    /* Define centers of refinement and point creation. */
    g->a[0] = 0.2;
    g->a[1] = 0.4;
    g->a[2] = 0.4;
    g->b[0] = 0.7;
    g->b[1] = 0.55;
    g->b[2] = 0.55;
    g->c[0] = 0.3;
    g->c[1] = 0.8;
    g->c[2] = 0.8;

    /*
     * Run example.
     */
    run (g);
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
