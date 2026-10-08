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

/** \file t8_mesh_step4_partition_balance_ghost.cxx
 * This is step4 of the t8code mesh handle tutorials.
 * Therefore, this is the same as general/t8_step4_partition_balance_ghost.cxx but using the mesh handle interface instead of the forest 
 * interface.
 * After generating a coarse mesh (step1), building a uniform mesh
 * on it (step2) and adapting this mesh (step3)
 * we will now learn how to partition and balance a mesh and how to generate a layer of ghost elements. 
 * Note: Normally, you would call set_adapt(), set_partition(), set_balance() and set_ghost() and then commit() only once, 
 * as shown in step 5. commit() applies the operations in the correct order automatically, regardless of the order of the set_*() calls. 
 * In this tutorial, we commit after each operation instead, so that we can see what each of them does and 
 * compare the resulting vtk files one after another.
 */

#include <t8.h>                                 /** General t8code header. Always include this. */
#include <mesh_handle/mesh.hxx>                 /** General mesh header. Always needed for mesh_handle code. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic cmesh to mesh_handle conversions. */
#include <mesh_handle/concepts.hxx> /** Include this to use c++ concepts related to the mesh handle. This can be used to constrain the template parameters to only allow mesh handle classes. */
#include <t8_types/t8_vec.hxx>      /** t8 vector dataclass. */
#include "t8_mesh_tutorials_common.hxx" /** Adaption function definition used for this tutorial. */
#include <memory>

using mesh_type = t8_mesh_handle::
  mesh<>; /**< Mesh class used in this tutorial. We define it globally to get rid of the function templates to simplify the code. */

/** Helper function to print the total number of elements in the mesh after each step.
 *  \param [in] mesh  The mesh handle to get the number of elements from.
 *  \param [in] stage The stage of the mesh (e.g. "Initial mesh", "Adapted mesh", etc.) to print in the output.
 *  \param [in] prefix The prefix for the filename of the exported vtk files. 
*/
void
print_stats_and_export (mesh_type& mesh, const char* stage, const char* prefix)
{
  t8_locidx_t local_elements = mesh.get_num_local_elements ();
  t8_gloidx_t global_elements = mesh.get_num_global_elements ();

  t8_global_productionf (" [mesh_step4] === %s === \n", stage);
  t8_global_productionf (" [mesh_step4] Local elements on the root process: %i \n", local_elements);
  t8_global_productionf (" [mesh_step4] Total elements: %li \n", global_elements);

  /* Writing the mesh to vtu and pvtu files, using the extended version of the function to ensure additional data like ghost elements, treeid etc. to be written into the files. */
  t8_mesh_handle::write_mesh_to_vtk_ext (mesh, prefix, 0, nullptr, true, true, true, true, true, false, false);
}

/** Helper function to adapt a given mesh using the predefined adaption callback function.
 *  \param [in] mesh  The initial mesh to adapt.
 *  \param [in] adapt_params  The adaptation parameters to use for the adaptation.
*/

void
step4_adapt_mesh (mesh_type& mesh, const adapt_data& adapt_params)
{
  /* Setting the adapt-flag with our adapt_callback_sphere function from step 3 and the adapt_params. Both can be found in the file \ref t8_mesh_tutorials_common.hxx. */
  mesh.set_adapt (mesh_type::mesh_adapt_callback_wrapper<adapt_data> (&adapt_callback_sphere<mesh_type>, adapt_params));
  /* Committing the mesh. */
  mesh.commit ();
}

/** Helper function to partition and balance a given mesh.
 *  \param [in] mesh  The initial mesh to partition and balance.
*/

void
step4_partition_balance_mesh (mesh_type& mesh)
{
  /* Setting partition flag.*/
  mesh.set_partition ();

  /* Setting balancing flag. */
  mesh.set_balance ();

  /* Committing the mesh. */
  mesh.commit ();
}

/** Helper function to add a layer of ghost elements to an initial mesh.
 *  \param [in] mesh  The initial mesh to add the ghost elements to.
*/
void
step4_ghost_mesh (mesh_type& mesh)
{
  /* Set flag such that ghost layer is created on commit. */
  mesh.set_ghost ();

  /* Committing the mesh. */
  mesh.commit ();
}

/** Entry point of the program. */
int
main (int argc, char** argv)
{
  /* Setting the file prefixes for the vtk files. */
  const char* prefix_initial = "step4_initial_mesh";
  const char* prefix_adapt = "step4_adapted_mesh";
  const char* prefix_partition_balance = "step4_partitioned_balanced_mesh";
  const char* prefix_ghost = "step4_ghost_mesh";
  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Print a starting message. */
  t8_global_productionf (" [mesh_step4] \n");
  t8_global_productionf (" [mesh_step4] Hello, this is step 4 of the t8code mesh_handle tutorials.\n");
  t8_global_productionf (" [mesh_step4] In this example we will create a mesh, adapt, partition, balance "
                         "and create a ghost layer on it. \n");
  t8_global_productionf (" [mesh_step4] \n");

  /* The initial uniform refinement level. */
  const int uniform_level = 3;

  /* Parameters for the adaption step. */
  adapt_data adapt_params = { { 0.5, 0.5, 1.0 }, 0.2, 0.4 };

  /**
   * INITIAL MESH
  */

  t8_global_productionf (" [mesh_step4] \n");
  t8_global_productionf (" [mesh_step4] Creating initial mesh.\n");
  t8_global_productionf (" [mesh_step4] \n");
  { /** Mesh scope begin. */
    /* Creating the initial mesh with uniform refinement. */
    auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_type> (uniform_level, comm);

    /* Printing the mesh information. */
    print_stats_and_export (*mesh, "Initial mesh", prefix_initial);

    /** 
     * ADAPT MESH
     */

    t8_global_productionf (" [mesh_step4] \n");
    t8_global_productionf (" [mesh_step4] Adapt the mesh.\n");
    t8_global_productionf (" [mesh_step4] \n");

    /** Call adaption helper function. */
    step4_adapt_mesh (*mesh, adapt_params);

    /* Printing the mesh information. */
    print_stats_and_export (*mesh, "Adapted mesh", prefix_adapt);

    /**
     * PARTITION, BALANCE MESH
     */

    t8_global_productionf (" [mesh_step4] \n");
    t8_global_productionf (" [mesh_step4] Partition and balance the mesh.\n");
    t8_global_productionf (" [mesh_step4] \n");

    /** Adapting the mesh from above a second time to see a difference when balancing. */
    step4_adapt_mesh (*mesh, adapt_params);

    /** Call partitioning and balancing helper function. */
    step4_partition_balance_mesh (*mesh);

    /* Printing the mesh information. */
    print_stats_and_export (*mesh, "Partitioned and Balanced mesh", prefix_partition_balance);

    /**
     * GHOST LAYER
     */

    t8_global_productionf (" [mesh_step4] \n");
    t8_global_productionf (" [mesh_step4] Creating ghost layer for mesh.\n");
    t8_global_productionf (" [mesh_step4] \n");

    /** Call ghost helper function. */
    step4_ghost_mesh (*mesh);

    /* Printing the mesh information. */
    print_stats_and_export (*mesh, "Ghost mesh", prefix_ghost);
    int ghost_elements = mesh->get_num_ghosts ();
    t8_global_productionf (" [mesh_step4] Number of ghost elements: %i \n", ghost_elements);

    t8_global_productionf (" [mesh_step4] \n");
    t8_global_productionf (" [mesh_step4] Finished all steps successfully.\n");
    t8_global_productionf (" [mesh_step4] \n");
  } /** Mesh scope end. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
