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

/** \file t8_mesh_element_data.cxx
 * This is step4 of the t8code mesh handle tutorials.
 * Therefore, this is the same as general/t8_step4_partition_balance_ghost.cxx but using the mesh handle interface instead of the forest 
 * interface.
 * After generating a coarse mesh (step1), building a uniform mesh
 * on it (step2) and adapting this mesh (step3)
 * we will now learn how to control the mesh creation in more detail,
 * how to partition and balance a mesh and how to generate a layer of ghost elements.
 */

#include <t8.h>                                 /** General t8code header. Always include this. */
#include <mesh_handle/mesh.hxx>                 /** General mesh header. Always needed for mesh_handle code. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic cmesh to mesh_handle conversions. */
#include <mesh_handle/concepts.hxx> /** Include this to use c++ concepts related to the mesh handle. This can be used to constraint the template parameters to only allow mesh handle classes. */
#include <t8_types/t8_vec.hxx>      /** t8 vector dataclass. */
#include "t8_mesh_tutorials_common.hxx" /** Default adaption function. */
#include <memory>
#include <iostream>

using mesh_type = t8_mesh_handle::
  mesh<>; /**< Mesh class used in this tutorial. We define it globally to get rid of the function templates to simplify the code. */

/** Helper function to print the total number of elements in the mesh after each step.
 *  \param mesh  The mesh handle to get the number of elements from.
 *  \param stage The stage of the mesh (e.g. "Initial mesh", "Adapted mesh", etc.) to print in the output.
 *  \param comm  The MPI communicator to use for the reduction and printing.
*/
void
print_mesh_stats (const std::unique_ptr<mesh_type>& mesh, const char* stage, sc_MPI_Comm comm)
{
  int local_elements = mesh->get_num_local_elements ();
  int global_elements = 0;
  MPI_Allreduce (&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, comm);

  t8_global_productionf (" [mesh_step4] === %s === \n", stage);
  t8_global_productionf (" [mesh_step4] Total elements: %i \n", global_elements);
}

/** Helper function to create an adapted mesh from an initial mesh.
 *  \param mesh  The initial mesh to adapt.
 *  \param adapt_params  The adaptation parameters to use for the adaptation.
*/
void
create_adapted_mesh (std::unique_ptr<mesh_type>& mesh, const adapt_data& adapt_params)
{
  /* Adapting the mesh once with our adapt_callback_sphere function from step 3 and the adapt_params. Both can be found in the file \ref t8_mesh_tutorials_common.hxx. */
  mesh->set_adapt (
    mesh_type::template mesh_adapt_callback_wrapper<adapt_data> (&adapt_callback_sphere<mesh_type>, adapt_params));
  /* Committing the adapted mesh. */
  mesh->commit ();
}

/** Helper function to create a partitioned and balanced mesh from an initial mesh.
 *  \param mesh  The initial mesh to adapt.
*/
void
create_partitioned_balanced_mesh (const std::unique_ptr<mesh_type>& mesh)
{
  /* Calculate partition information.*/
  mesh->set_partition ();

  /* Calculate balancing information. */
  mesh->set_balance ();

  /* Committing the mesh --> modifying the mesh according to the precalculated information. */
  mesh->commit ();
}

/** Helper function to create a mesh with ghosts from an initial mesh.
 *  \param mesh  The initial mesh to adapt.
*/
void
create_ghost_mesh (const std::unique_ptr<mesh_type>& mesh)
{
  /* Creating the ghost layers. */
  mesh->set_ghost ();

  /* Committing the ghost mesh. */
  mesh->commit ();
}

/** Entry point of the program. */
int
main (int argc, char** argv)
{

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
  t8_global_productionf (" [t8 step 4 Mesh handle] \n");
  t8_global_productionf (
    " [t8 step 4 Mesh handle] Hello, this is the mesh adaptation example of t8code using the mesh handle.\n");
  t8_global_productionf (" [t8 step 4 Mesh handle] In this example we will create a mesh, adapt, partition, balance "
                         "and create a ghost layer on it. \n");
  t8_global_productionf (" [t8 step 4 Mesh handle] \n");

  /* The initial uniform refinement level. */
  int uniform_level = 3;

  /* Parameters for the adaption step. */
  struct adapt_data adapt_params = { { 0.5, 0.5, 1.0 }, 0.2, 0.4 };

  /**
   * INITIAL MESH
  */

  t8_global_productionf (" [t8 step 4 Mesh handle] \n");
  t8_global_productionf (" [t8 step 4 Mesh handle] Creating initial mesh.\n");
  t8_global_productionf (" [t8 step 4 Mesh handle] \n");
  { /** Mesh scope begin. */
    /* Creating the initial mesh with uniform refinement. */
    auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_type> (uniform_level, comm);

    /* Printing the mesh information. */
    print_mesh_stats (mesh, "Initial mesh", comm);

    /* Writing the mesh to vtu and pvtu files, using the extended version of the function to ensure additional data like ghost elements, treeid etc. to be written into the files. */
    t8_mesh_handle::write_mesh_to_vtk_ext (*mesh, "initial_mesh.vtu", 0, nullptr, true, true, true, true, true, false,
                                           false);

    /** 
   * ADAPTED MESH
  */

    t8_global_productionf (" [t8 step 4 Mesh handle] \n");
    t8_global_productionf (" [t8 step 4 Mesh handle] Creating adapted mesh.\n");
    t8_global_productionf (" [t8 step 4 Mesh handle] \n");

    /** Call adaption helper function. */
    create_adapted_mesh (mesh, adapt_params);

    /* Printing the mesh information. */
    print_mesh_stats (mesh, "Adapted mesh", comm);

    /* Writing the mesh to vtu and pvtu files using the extended version of the function. */
    t8_mesh_handle::write_mesh_to_vtk_ext (*mesh, "adapted_mesh.vtu", 0, nullptr, true, true, true, true, true, false,
                                           false);

    /**
   * PARTITIONED, BALANCED MESH
  */

    t8_global_productionf (" [t8 step 4 Mesh handle] \n");
    t8_global_productionf (" [t8 step 4 Mesh handle] Creating partitioned and balanced mesh.\n");
    t8_global_productionf (" [t8 step 4 Mesh handle] \n");

    /** Adapting the mesh from above a second time to see a difference when balancing. */
    create_adapted_mesh (mesh, adapt_params);

    /** Call partitioning and balancing helper function. */
    create_partitioned_balanced_mesh (mesh);

    /* Printing the mesh information. */
    print_mesh_stats (mesh, "Partitioned and Balanced mesh", comm);

    /* Writing the mesh to vtu and pvtu files using the extended version of the function. */
    t8_mesh_handle::write_mesh_to_vtk_ext (*mesh, "partition_balance_mesh.vtu", 0, nullptr, true, true, true, true,
                                           true, false, false);

    /**
   * GHOST MESH
  */

    t8_global_productionf (" [t8 step 4 Mesh handle] \n");
    t8_global_productionf (" [t8 step 4 Mesh handle] Creating ghost layer for mesh.\n");
    t8_global_productionf (" [t8 step 4 Mesh handle] \n");

    /** Call ghost helper function. */
    create_ghost_mesh (mesh);

    /* Printing the mesh information. */
    print_mesh_stats (mesh, "Ghost mesh", comm);

    /* Writing the mesh to vtu and pvtu files using the extended version of the function. */
    t8_mesh_handle::write_mesh_to_vtk_ext (*mesh, "ghost_mesh.vtu", 0, nullptr, true, true, true, true, true, false,
                                           false);

    t8_global_productionf (" [t8 step 4 Mesh handle] \n");
    t8_global_productionf (" [t8 step 4 Mesh handle] Finished all steps successfully.\n");
    t8_global_productionf (" [t8 step 4 Mesh handle] \n");
  } /** Mesh scope end. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
