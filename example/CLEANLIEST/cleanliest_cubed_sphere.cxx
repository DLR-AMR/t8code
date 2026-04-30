/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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
#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_functions.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <t8.h>                                 /* General t8code header, always include this. */
#include <t8_cmesh/t8_cmesh.h>                  /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>         /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>        /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>             /* save forest */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_types/t8_vec.hxx>                  /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_iterate.h>        /* For the search algorithm. */
#include <tutorials/general/t8_step3.h>         /* Example forest adaptation from step 3 */
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_profiling.h>

#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>

#include "t8_adapt_to_point_sources.h"
#include "t8_point_source_parser.h"

#include <cmath>
#include <numbers>

/** Create an array of a given number of particles on the root process
 * and broadcast it to all other processes.
 * \param [in] num_particles  The number of particles to create.
 * \param [in] seed           The seed to be used for the random number generator.
 * \param [in] comm           MPI communicator to specify on which processes we create this array.
 */
static sc_array *
t8_spherical_shell_random_particles (size_t num_particles, unsigned int seed, double inner_radius, 
                                     double shell_thickness, sc_MPI_Comm comm)
{
  /* Specify lower and upper bounds for the coordinates in each dimension. */
  double boundary_low[3] = { 0.0, 0.0, 0.0 };
  double boundary_high[3] = { 1.0, 1.0, 1.0 };
  int mpirank;
  int mpiret;
  sc_array *particles;
  const unsigned int nsd = 3;

  /* Get the MPI rank. */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Create an array for num_particles many particles. */
  particles = sc_array_new_count (sizeof (t8_point_source_t), num_particles);

  /* We build the array on rank 0 and broadcast it to the other ranks.
   * This ensures that all ranks have the same randomly generated particles. */
  if (mpirank == 0) {
    /* Rank 0 fills this array with random particles. */
    size_t iparticle;
    srand (seed);

    double rand_number[3] = {};

    for (iparticle = 0; iparticle < num_particles; ++iparticle) {
      /* Get this particle's pointer. */
      t8_point_source_t *particle = (t8_point_source_t *) sc_array_index_int (particles, iparticle);

      // double radius = 0.0;
      for (unsigned int dim = 0; dim < nsd; ++dim) {
        rand_number[dim] = (double) rand () / RAND_MAX;
        // radius += rand_number[dim] * rand_number[dim];
      }
      // radius = sqrt(radius);

      // Transform random number into (some pseudo-random) spherical coordinates 
      // (1) Distance r
      double r = inner_radius + rand_number[0] * shell_thickness;
      // (2,3) Angles: theta and phi
      double theta = rand_number[1] * std::numbers::pi;
      double phi   = rand_number[2] * std::numbers::pi * 2.0;

      // Convert them to Cartesian coordinates
      particle->coordinates[0] = r * sin(theta) * cos(phi);
      particle->coordinates[1] = r * sin(theta) * sin(phi);
      particle->coordinates[2] = r * cos(theta);

      /* Initialize the is_inside_partition flag. */
      particle->is_inside_partition = 0;

      t8_debugf ("Coordinates: %.3f %.3f %.3f\n", particle->coordinates[0], particle->coordinates[1],
                 particle->coordinates[2]);
    }
  }
  /* Broadcast this array to all other processes. */
  mpiret = sc_MPI_Bcast (particles->array, sizeof (t8_point_source_t) * num_particles, sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);

  t8_global_productionf (
    " [search] Created %zd random particles inside the box [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f].\n", num_particles,
    boundary_low[0], boundary_high[0], boundary_low[1], boundary_high[1], boundary_low[2], boundary_high[2]);

  return particles;
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  size_t num_particles;
  int level;
  int num_layers;
  int num_horizontal_trees_per_quarter;
  const unsigned seed = 0;

  sc_options_t *opt;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is some cubed sphere test for CLEANLIEST.\n");
  opt = sc_options_new (argv[1]);
  sc_options_add_int (opt, 'l', "level", &level, 5, "The level of the forest.");
  sc_options_add_size_t (opt, 'p', "num_particles", &num_particles, 20, "The number of particles.");
  sc_options_add_int (opt, 'v', "vertical_layers", &num_layers, 5, "The number of shell layers (vertical direction) of the cmesh.");
  sc_options_add_int (opt, 'h', "horizontal_trees", &num_horizontal_trees_per_quarter, 5, "The number of horizontal trees per quarter of the spherical shell.");
  sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT, opt, argc, argv);

  t8_cmesh_t cmesh;
  t8_forest_t forest;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);
  comm = sc_MPI_COMM_WORLD;

  double earth_radius = 6.0e+3;
  double atmosphere_height = 50.0;


  // Create cmesh: Spherical shell representing the atmosphere.
  cmesh = t8_cmesh_new_cubed_spherical_shell(earth_radius, atmosphere_height, num_horizontal_trees_per_quarter, num_layers, comm);

  // Create uniform forest
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, comm);

  // Write uniform forest to vtk
  t8_forest_write_vtk (forest, "t8_cubed_sphere_uniform");

  // Parse point sources from CSV file into std::vector
  auto parsedPointSources = t8_parse_point_sources_from_file("global_emission_sources.csv");

  // Convert std::vector<PointSourceH2> to sc_array of t8_point_source_t:
  sc_array *point_sources;
  int num_point_sources = parsedPointSources.size();
  point_sources = sc_array_new_count (sizeof (t8_point_source_t), num_point_sources);
  for( unsigned int isrc=0; isrc < parsedPointSources.size(); isrc++)
  {
    // Access SC array entry as pointer
    t8_point_source_t *cur_point_source = (t8_point_source_t *) sc_array_index_int (point_sources, isrc);

    // Copy coordinates
    cur_point_source->coordinates[0] = parsedPointSources[isrc].coordinates.at(0);
    cur_point_source->coordinates[1] = parsedPointSources[isrc].coordinates.at(1);
    cur_point_source->coordinates[2] = parsedPointSources[isrc].coordinates.at(2);

    // t8_glo

    /* Initialize the is_inside_partition flag. */
    cur_point_source->is_inside_partition = 0;
  }

  /* Broadcast this array to all other processes. */
  mpiret = sc_MPI_Bcast (point_sources->array, sizeof (t8_point_source_t) * num_point_sources, sc_MPI_BYTE, 0, comm);
  SC_CHECK_MPI (mpiret);


  // /* Create an array with random particles. */
  // sc_array_t *particles = t8_spherical_shell_random_particles (num_particles, seed, earth_radius, atmosphere_height, comm);

  // Adapt forest
  forest = t8_adapt_to_point_sources (forest, point_sources, "cubed_spherical_shell_adapted", false);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);

  sc_array_destroy (point_sources);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
