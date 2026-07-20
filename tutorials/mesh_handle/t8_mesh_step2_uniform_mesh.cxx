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

/** \file t8_mesh_step2_uniform_mesh.cxx
 * This is step 2 of the t8code mesh handle tutorials.
 * Therefor, this is the same as general/t8_step2_uniform_forest.cxx but using the mesh handle interface instead of the forest
 * interface.
 * After we learned how to create a cmesh in step1, we will
 * now build our first partitioned mesh, get its local and global
 * element count, and output it into .vtu files.
 *
 * When we create a mesh from a coarse mesh using the mesh handle interface, the mesh will always be
 * uniform (every element has the same refinement level) and can then be adapted
 * later (see the following steps).
 * Together with the cmesh, we also need a refinement scheme. This scheme tells the
 * mesh how elements of each shape (t8_eclass_t) are refined, what their neighbor
 * are etc.
 * The default scheme in t8_schemes/t8_default/t8_default.hxx provides an implementation for
 * all element shapes that t8code supports (with pyramids currently under construction).
*/
#include <t8.h>                                 /** General t8code header, always include this. */
#include <mesh_handle/mesh.hxx>                 /** General Mesh Header, always needed for mesh_handle code. */
#include <t8_cmesh/t8_cmesh.h>                  /** cmesh definition and basic interface. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic Cmesh to mesh_handle conversions. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <t8_schemes/t8_default/t8_default.hxx> /** default refinement scheme. */
#include <string>

/** Builds cmesh of 2 prisms that build up a unit cube.
 * See step1 for a detailed description.
 * \param [in] comm   MPI Communicator to use.
 * \return            The coarse mesh.
 */
static t8_cmesh_t
t8_step2_build_prismcube_coarse_mesh (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;

  /* Build a coarse mesh of 2 prisms that form a cube. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_PRISM, comm, 0, 0, 0);
  t8_global_productionf (" [tutorial] Constructed coarse mesh with 2 prisms.\n");

  return cmesh;
}

/** Build a uniform mesh on a cmesh using the default refinement scheme. 
 * \param [in] comm   MPI Communicator to use.
 * \param [in] cmesh  The coarse mesh to build the uniform mesh on.
 * \param [in] level  The initial uniform refinement level.
 * \return            A uniform mesh with the given refinement level that is
 *                    partitioned across the processes in \a comm.
 */
static std::unique_ptr<t8_mesh_handle::mesh<t8_mesh_handle::all_cache_element_competences>>
t8_step2_build_uniform_mesh (sc_MPI_Comm comm, t8_cmesh_t cmesh, int level)
{
  const t8_scheme *scheme = t8_scheme_new_default (); /** Default refinement scheme. */

  /* Build the uniform mesh, it is automatically partitioned among the processes. */
  auto mesh = t8_mesh_handle::handle_new_uniform<t8_mesh_handle::mesh<t8_mesh_handle::all_cache_element_competences>> (
    cmesh, scheme, level, comm, false);

  t8_global_productionf (" [tutorial] Constructed uniform mesh with %d elements per tree.\n", 1 << (3 * level));

  return mesh;
}

int
main (int argc, char **argv)
{
  int mpiret;
  sc_MPI_Comm comm;

  /** File prefix for our vtk files. */
  const char *prefix = "t8_step2_uniform_mesh";
  /** Uniform refinement level of the mesh. */
  const int level = 3;
  t8_locidx_t local_num_elements;
  t8_gloidx_t global_num_elements;

  /** Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /** Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /** Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /** Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);

  /** Print a message on the root process. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Hello, this is the step2 example of t8code using the mesh handle.\n");
  t8_global_productionf (" [tutorial] In this example we build our first uniform mesh and output it to vtu files.\n");
  t8_global_productionf (" [tutorial] \n");

  /** We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /** Create the cmesh. */
  t8_cmesh_t cmesh = t8_step2_build_prismcube_coarse_mesh (comm);

  /** Build the uniform mesh. */
  auto mesh = t8_step2_build_uniform_mesh (comm, cmesh, level);
  /** Get the number of local elements. */
  local_num_elements = mesh->get_num_local_elements ();
  /** Get the number of global elements. */
  global_num_elements = mesh->get_num_global_elements ();

  /** Print information on the mesh. */
  t8_global_productionf (" [tutorial] Created uniform mesh.\n");
  t8_global_productionf (" [tutorial] Refinement level:\t\t\t%i\n", level);
  t8_global_productionf (" [tutorial] Local number of elements:\t\t%i\n", local_num_elements);
  t8_global_productionf (" [tutorial] Global number of elements:\t%" T8_GLOIDX_FORMAT "\n", global_num_elements);

  /** Write mesh to vtu files. */
  t8_mesh_handle::write_mesh_to_vtk (*mesh, prefix);
  t8_global_productionf (" [tutorial] Wrote mesh to vtu files:\t%s*\n", prefix);

  /** Destroy the mesh. */
  mesh.reset ();
  t8_global_productionf (" [tutorial] Destroyed mesh.\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
