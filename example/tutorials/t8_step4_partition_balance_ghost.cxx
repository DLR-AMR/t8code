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

/* See also: https://github.com/holke/t8code/wiki/Step-4---Partition,-Balance,-Ghost
 *
 * This is step4 of the t8code tutorials.
 * After generating a coarse mesh (step1), building a uniform forest
 * on it (step2) and adapting this forest (step3) 
 * we will now lear how to control the forest creation in more detail,
 * how to partition and balance a forest and how to generate a layer of ghost elements.
 * 
 * Partition: Each forest is distributed among the MPI processes. Partitioning a forest means
 *            to change this distribution in order to maintain a load balance. After we 
 *            applied partition to a forest the new forest will have equal numbers of elements
 *            on each process (+- 1).
 *            In our example we start with a uniform forest from a cmesh. This forest is constructed
 *            such that each process has the same number of elements. 
 *            We then adapt this forest, thus refining and coarsening its elements. This changes the
 *            number of elements on each process and will almost always result in a load imbalance.
 *            You can verify this yourself by printing the process local number on each process for the 
 *            adapted forest (for example with t8_productionf).
 *            In order to reestablish a load balance, we will construct a new forest from the adapted
 *            one via the partition feature.
 * 
 * Ghost:     Many applications require a layer of ghost elements. Ghost element of a process are 
 *            those elements that are not local to this process but have a (face) neighbor element that is.
 *            Telling a forest to create a ghost layer will gather all necessary information and give
 *            us access to the ghost elements.
 *            t8code does not create a layer of ghost elements by default, thus our initial uniform forest
 *            and the adapted forest will not have one.
 * 
 * Balance:   A forest fulfills the balance property if (and only if) for each element its (face) neighbors
 *            have refinement level at most +1 or -1 in comparison to the element's level.
 *            The balance property is often broken after adaptation. The Balance algorithm iterates through
 *            the mesh and restores the balance property by refining elements if necessary.
 *            Balance will never coarsen any elements.
 *            t8code does not balance a forest by default.
 *            Our initial uniform forest is balanced since every element has the same refinement level l.
 *            Our adaptation criterion is such that (usually) after one step of adaptation the forest will still be balanced,
 *            since we have level l+1 elements in the inner sphere, level l elements in the middle and level l+1 elements
 *            outside. This may not be the case for very small initial refinement levels levels or with different radius thresholds (why don't you try it out?).
 *            Therefore, in this example we will apply the adaptation two times, resulting in level l+2 elements
 *            in the inner sphere, level l elements in the middle and level l - 2 element in the outher sphere
 *            (and probably some, but not many, level l-1 and level l+1 elements in between).
 *            This forest will be unbalanced and we will then apply the balance routine to it.
 *            Note that balance changes the local number of elements and thus may also change the load balance
 *            and hence require repartitioning.
 *            Balance is usually the most expensive of t8code's mesh manipulation algorithms.
 * 
 * How you can experiment here:
 *    Partition:
 *         - Test the program with different numbers of processes and compare the local
 *           number of elements before and after partitioning.
 *         - Open the vtu files before and after partitioning and color them by mpirank in order
 *           to see how the distribution of elements changes.
 *    Ghost:
 *         - We actually print the ghost elements in the vtu files and you can look at them.
 *           Ghost elements have treeid = -1 in the vtu files. To view the ghosts of a single process (say process 1)
 *           you can set a threshhold on mpirank to view only the elements of process 1 and then set a threshold
 *           on treeid to only view those elements with treeid = -1.
 *    Balance:
 *         - View the unbalanced and balanced forest vtu files and compare them.
 *           You can color the elements by level.       
 *         - Exercise: Apply balance to forest_adapt and not to the twice adapted forest.
 *         - Set the no_repartition flag of t8_forest_set_balance to true and observe how the
 *           partition of the resulting forest changes.
 *  */

#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <t8_schemes/t8_default_cxx.hxx>        /* default refinement scheme. */
#include <example/tutorials/t8_step3.h>

T8_EXTERN_C_BEGIN ();

  /* So far we have seen t8_forest_new_* functions to create forests.
   * These directly returned a new forest.
   * However, t8code offers us more control over the creation of forests.
   * For example we can control wether or not a forest should have a ghost layer,
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

/* In this function we create a new forest that repartitions a given forest
 * and has a layer of ghost elements. 
 */
static t8_forest_t
t8_step4_partition_ghost (t8_forest_t forest)
{
  t8_forest_t         new_forest;

  /* Check that forest is a committed, that is valid and usable, forest. */
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Initialize */
  t8_forest_init (&new_forest);

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
  t8_forest_set_partition (new_forest, forest, 0);
  /* Tell the new_forest to create a ghost layer.
   * This will gather those face neighbor elements of process local element that reside
   * on a different process.
   * 
   * We currently support ghost mode T8_GHOST_FACES that creates face neighbor ghost elements
   * and will in future also support other modes for edge/vertex neighbor ghost elements.
   */
  t8_forest_set_ghost (new_forest, 1, T8_GHOST_FACES);
  /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
  t8_forest_commit (new_forest);

  return new_forest;
}

/* In this function we adapt a forest as in step3 and balance it. 
 * In our main program the input forest is already adapted and then the resulting twice adapted forest will be unbalanced.
 */
static t8_forest_t
t8_step4_balance (t8_forest_t forest)
{
  t8_forest_t         balanced_forest;
  /* Adapt the input forest. */
  t8_forest_t         unbalanced_forest = t8_step3_adapt_forest (forest);

  /* Output to vtk. */
  t8_forest_write_vtk (unbalanced_forest, "t8_step4_unbalanced_forest");
  t8_global_productionf
    (" [step4] Wrote unbalanced forest to vtu files: %s*\n",
     "t8_step4_unbalanced_forest");

  /* Initialize new forest. */
  t8_forest_init (&balanced_forest);
  /* Specify that this forest should result from balancing unbalanced_forest.
   * The last argument is the flag 'no_repartition'.
   * Since balancing will refine elements, the load-balance will be broken afterwards.
   * Setting this flag to false (no_repartition = false -> yes repartition) will repartition
   * the forest after balance, such that every process has the same number of elements afterwards.
   */
  t8_forest_set_balance (balanced_forest, unbalanced_forest, 0);
  /* Commit the forest. */
  t8_forest_commit (balanced_forest);

  return balanced_forest;
}

int
t8_step4_main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  /* The prefix for our output files. */
  const char         *prefix_uniform = "t8_step4_uniform_forest";
  const char         *prefix_adapt = "t8_step4_adapted_forest";
  const char         *prefix_partition_ghost =
    "t8_step4_partitioned_ghost_forest";
  const char         *prefix_balance = "t8_step4_balanced_forest";
  /* The uniform refinement level of the forest. */
  const int           level = 3;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* Print a message on the root process. */
  t8_global_productionf (" [step4] \n");
  t8_global_productionf
    (" [step4] Hello, this is the step4 example of t8code.\n");
  t8_global_productionf
    (" [step4] In this example we will explain the forest creation in more details.\n");
  t8_global_productionf
    (" [step4] We will partition a forest and create a ghost layer and we will also balance a forest.\n");
  t8_global_productionf (" [step4] \n");

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  /*
   * Setup.
   * Build cmesh and uniform forest.
   */

  t8_global_productionf (" [step4] \n");
  t8_global_productionf
    (" [step4] Creating an adapted forest as in step3.\n");
  t8_global_productionf (" [step4] \n");
  /* Build a cube cmesh with tet, hex, and prism trees. */
  cmesh = t8_cmesh_new_hypercube_hybrid (comm, 0, 0);
  t8_global_productionf (" [step4] Created coarse mesh.\n");
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);

  /* Print information of the forest. */
  t8_step3_print_forest_information (forest);

  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_uniform);
  t8_global_productionf
    (" [step4] Wrote uniform level %i forest to vtu files: %s*\n", level,
     prefix_uniform);

  /*
   *  Adapt the forest.
   */

  /* Adapt the forest. We can reuse the forest variable, since the new adapted
   * forest will take ownership of the old forest and destroy it.
   * Note that the adapted forest is a new forest, though. */
  forest = t8_step3_adapt_forest (forest);

  /* Print information of our new forest. */
  t8_step3_print_forest_information (forest);

  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_adapt);
  t8_global_productionf (" [step4] Wrote adapted forest to vtu files: %s*\n",
                         prefix_adapt);

  /*
   * Partition and create ghost elements.
   */

  t8_global_productionf (" [step4] \n");
  t8_global_productionf
    (" [step4] Repartitioning this forest and creating a ghost layer.\n");
  t8_global_productionf (" [step4] \n");
  forest = t8_step4_partition_ghost (forest);
  t8_global_productionf
    (" [step4] Repartitioned forest and built ghost layer.\n");
  t8_step3_print_forest_information (forest);
  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_partition_ghost);

  /*
   * Balance
   */

  t8_global_productionf (" [step4] \n");
  t8_global_productionf
    (" [step4] Creating an unbalanced forest and balancing it.\n");
  t8_global_productionf (" [step4] \n");
  forest = t8_step4_balance (forest);
  t8_global_productionf (" [step4] Balanced forest.\n");
  t8_step3_print_forest_information (forest);
  /* Write forest to vtu files. */
  t8_forest_write_vtk (forest, prefix_balance);
  t8_global_productionf (" [step4] Wrote balanced forest to vtu files: %s*\n",
                         prefix_balance);

  /*
   * clean-up
   */

  /* Destroy the forest. */
  t8_forest_unref (&forest);
  t8_global_productionf (" [step4] Destroyed forest.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
