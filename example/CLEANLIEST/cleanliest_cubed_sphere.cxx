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

/** Create an array of a given number of particles on the root process
 * and broadcast it to all other processes.
 * \param [in] num_particles  The number of particles to create.
 * \param [in] seed           The seed to be used for the random number generator.
 * \param [in] comm           MPI communicator to specify on which processes we create this array.
 */
static sc_array *
t8_cubed_sphere_random_particles (size_t num_particles, unsigned int seed, unsigned int nsd, sc_MPI_Comm comm)
{
  /* Specify lower and upper bounds for the coordinates in each dimension. */
  double boundary_low[3] = { 0.0, 0.0, 0.0 };
  double boundary_high[3] = { 1.0, 1.0, 1.0 };
  int mpirank;
  int mpiret;
  sc_array *particles;

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

      for (unsigned int dim = 0; dim < nsd; ++dim) {
        /* Create a random value between boundary_low[dim] and boundary_high[dim] */
        particle->coordinates[dim] = rand_number[dim] / sqrt (3);
        /* Initialize the is_inside_partition flag. */
        particle->is_inside_partition = 0;
      }

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
  const unsigned seed = 0;

  sc_options_t *opt;

  /* Print a message on the root process. */
  t8_global_productionf (" [search] \n");
  t8_global_productionf (" [search] Hello, this is some cubed sphere test for CLEANLIEST.\n");
  opt = sc_options_new (argv[1]);
  sc_options_add_int (opt, 'l', "level", &level, 5, "The level of the forest.");
  sc_options_add_size_t (opt, 'n', "num_particles", &num_particles, 20, "The number of particles.");
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

  // Create cmesh for cubed sphere
  cmesh = t8_cmesh_new_cubed_sphere (1.0, comm);

  // Create uniform forest
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, comm);

  // Write uniform forest to vtk
  t8_forest_write_vtk (forest, "t8_cubed_sphere_uniform");

  /* Create an array with random particles. */
  unsigned int nsd = t8_forest_get_dimension (forest);
  sc_array_t *particles = t8_cubed_sphere_random_particles (num_particles, seed, nsd, comm);

  // Adapt forest
  forest = t8_adapt_to_point_sources (forest, particles, "cubed_sphere_adapted");

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);

  sc_array_destroy (particles);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
