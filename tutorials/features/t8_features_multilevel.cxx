/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2024 the developers

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

/* See also: https://github.com/holke/t8code/wiki/Step-4---Partition,-Balance,-Ghost
 *
 * This is step4 of the t8code tutorials.
 * After generating a coarse mesh (step1), building a uniform forest
 * on it (step2) and adapting this forest (step3) 
 * we will now lear how to control the forest creation in more detail,
 * how to partition and balance a forest and how to generate a layer of ghost elements.
 * 
 */

#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>             /* A collection of exemplary cmeshes */
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                 /* save forest */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <tutorials/general/t8_step3.h>

/* So far we have seen t8_forest_new_* functions to create forests.
   * These directly returned a new forest.
   * However, t8code offers us more control over the creation of forests.
   * For example we can control whether or not a forest should have a ghost layer,
   * be balanced/partitioned from another forest, etc.
   * 
   * Usually, there are three steps involved in creating a forest:
   * 1. Initialize the forest with t8_forest_init.
   *      This function will prepare a forest by setting default values for its members and 
   *      initializing necessary structures.
   * 2. Set all properties that the forest should have.
   *      Here we can for example set a cmesh and refinement scheme or specify that the
   *      forest should be adapted/partitioned/balanced from another forest etc.
   *      See the t8_forest_set functions in t8_forest.h
   *      The order in which you call the set functions does not matter.
   * 3. Commit the forest with t8_forest_commit. 
   *      In this step the forest is actually created after setting all
   *      desired properties. The forest cannot be changed after it was committed.   
   * 
   * The t8_forest_new functions are just wrappers around this process.
   */

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  /* The prefix for our output files. */
  const char *prefix_uniform = "t8_multilevel_uniform_forest";
  const char *prefix_partition = "t8_multilevel_partitioned_forest";
  /* The uniform refinement level of the forest. */
  const int level = 3;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);

  /* Print a message on the root process. */
  t8_global_productionf (" [feature multilevel] \n");
  t8_global_productionf (" [feature multilevel] Hello, this is the multilevel example of t8code.\n");
  t8_global_productionf (" [feature multilevel] In this example we will explain the multileel forest creation in more detail.\n");
  t8_global_productionf (
    " [feature multilevel] We will create a uniform multilevel forest and partition it.\n");
  t8_global_productionf (" [feature multilevel] \n");

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /*
   * Setup.
   * Build cmesh and uniform forest.
   */

  t8_global_productionf (" [feature multilevel] \n");
  t8_global_productionf (" [feature multilevel] Creating a uniform multilevel forest.\n");
  t8_global_productionf (" [feature multilevel] \n");
  /* Build a cube cmesh with one hex. */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
  t8_global_productionf (" [feature multilevel] Created coarse mesh.\n");
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, 1, comm);

  /* Print information of the forest. */
  //t8_step3_print_forest_information (forest);

  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_uniform);
  t8_global_productionf (" [feature multilevel] Wrote uniform level %i forest to vtu files: %s*\n", level, prefix_uniform);

  
  /*
   * Partition the forest.
   */

  t8_global_productionf (" [feature multilevel] \n");
  t8_global_productionf (" [feature multilevel] Repartitioning the forest .\n");
  t8_global_productionf (" [feature multilevel] \n");
  t8_forest_t new_forest;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Initialize */
  //t8_forest_init (&new_forest);

  /* Tell the new_forest that is should partition the existing forest.
   * This will change the distribution of the forest elements among the processes
   * in such a way that afterwards each process has the same number of elements
   * (+- 1 if the number of elements is not divisible by the number of processes).
   * 
   * The third 0 argument is the flag 'partition_for_coarsening' which is currently not
   * implemented. Once it is, this will ensure that a family of elements will not be split
   * across multiple processes and thus one level coarsening is always possible (see also the
   * comments on coarsening in t8_step3).
   */
  //t8_forest_set_partition (new_forest, forest, 0);
  /* Tell the new_forest to create a ghost layer.
   * This will gather those face neighbor elements of process local element that reside
   * on a different process.
   * 
   * We currently support ghost mode T8_GHOST_FACES that creates face neighbor ghost elements
   * and will in future also support other modes for edge/vertex neighbor ghost elements.
   */
  //t8_forest_commit (new_forest);
  //t8_global_productionf (" [feature multilevel] Repartitioned forest.\n");
  //t8_step3_print_forest_information (forest);
  /* Write forest to vtu files. */
  //t8_forest_write_vtk_ext (forest, prefix_partition, 1, 1, 1, 1, 0, 0, 1, 0, NULL);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" [feature multilevel] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

