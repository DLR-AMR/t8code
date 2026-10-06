/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_extruded_hex.cxx
 * Example for the extruded hex scheme: A brick of hexes is refined around a vertical cylinder.
 * The hexes are only refined in x- and y-direction and always span the whole height of their tree, so the vertical
 * resolution is given by the number of trees stacked in z-direction.
 * The forest is balanced, partitioned, gets a ghost layer and is written to vtk.
 */

#include <t8.h>
#include <sc_options.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_schemes/t8_extruded/t8_extruded.hxx>
#include <t8_types/t8_vec.hxx>
#include <cmath>

/** The parameters of the refinement. */
struct t8_extruded_hex_adapt_data
{
  t8_2D_vec center; /**< The center of the cylinder in the x-y-plane. */
  double radius;    /**< The radius of the cylinder. */
  double width;     /**< Elements closer than this to the cylinder surface are refined. */
  int maxlevel;     /**< The maximum refinement level. */
};

/** Refine elements close to the surface of a vertical cylinder. */
static int
t8_extruded_hex_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                       const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                       [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                       t8_element_t *elements[])
{
  const t8_extruded_hex_adapt_data *data = (const t8_extruded_hex_adapt_data *) t8_forest_get_user_data (forest);
  if (scheme->element_get_level (tree_class, elements[0]) >= data->maxlevel) {
    return 0;
  }
  t8_3D_vec centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data ());
  t8_2D_vec centroid2D ({ centroid[0], centroid[1] });
  const double distance = t8_dist (centroid2D, data->center);
  return std::abs (distance - data->radius) < data->width ? 1 : 0;
}

int
main (int argc, char **argv)
{
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
#if T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_ESSENTIAL);
#endif

  int num_x, num_y, num_z, level, maxlevel, helpme;

  /* initialize command line argument parser */
  sc_options_t *opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'x', "num-x", &num_x, 4, "Number of trees in x-direction.");
  sc_options_add_int (opt, 'y', "num-y", &num_y, 4, "Number of trees in y-direction.");
  sc_options_add_int (opt, 'z', "num-z", &num_z, 3, "Number of trees in z-direction (number of layers).");
  sc_options_add_int (opt, 'l', "level", &level, 1, "Initial uniform refinement level.");
  sc_options_add_int (opt, 'm', "maxlevel", &maxlevel, 5, "Maximum refinement level.");

  const int parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme || parsed < 0 || num_x < 1 || num_y < 1 || num_z < 1 || level < 0 || maxlevel < level) {
    /* display help message and usage */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

    /* A brick of hexes, the z-direction is the extrusion direction. */
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_brick_3d (cmesh, num_x, num_y, num_z, 0, 0, 0, comm);
    SC_CHECK_ABORT (t8_cmesh_is_extrusion_compatible (cmesh), "The cmesh is not compatible with extruded hexes.");

    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_extruded (), level, 0, comm);
    t8_global_productionf ("Uniform forest of level %i with %lli elements.\n", level,
                           (long long) t8_forest_get_global_num_leaf_elements (forest));

    /* Refine around a vertical cylinder */
    t8_extruded_hex_adapt_data adapt_data;
    adapt_data.center[0] = num_x / 2.0;
    adapt_data.center[1] = num_y / 2.0;
    adapt_data.radius = 0.3 * std::min (num_x, num_y);
    adapt_data.width = 0.5;
    adapt_data.maxlevel = maxlevel;

    /* Set adapt data, partition, balance, and create a ghost layer. */
    t8_forest_t forest_adapt;
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, &adapt_data);
    t8_forest_set_adapt (forest_adapt, forest, t8_extruded_hex_adapt, 1);
    t8_forest_set_partition (forest_adapt, NULL, 0);
    t8_forest_set_balance (forest_adapt, NULL, 0);
    t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
    t8_forest_commit (forest_adapt);
    t8_global_productionf ("Adapted forest with %lli elements.\n",
                           (long long) t8_forest_get_global_num_leaf_elements (forest_adapt));

    t8_forest_write_vtk (forest_adapt, "t8_extruded_hex");
    t8_global_productionf ("Wrote the forest to t8_extruded_hex.\n");
    t8_forest_unref (&forest_adapt);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
