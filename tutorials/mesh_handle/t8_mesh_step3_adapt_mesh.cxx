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

/** \file t8_mesh_step3_adapt_mesh.cxx
 * This is step3 of the t8code mesh handle tutorials.
 * Therefore, this is the same as general/t8_step3_adapt_forest.cxx but using the mesh handle interface instead of the forest
 * interface.
 * After generating a coarse mesh (step1) and building a uniform mesh
 * on it (step2), we will now adapt (= refine and coarsen) the mesh
 * according to our own criterion.
 *
 * The geometry (coarse mesh) is again a cube, this time modelled with
 * 6 tetrahedra, 6 prisms and 4 cubes.
 * We refine an element if its midpoint is within a sphere of given radius
 * around the point (0.5, 0.5, 1) and we coarsen outside of a given radius.
 * We will use non-recursive refinement, that means that the refinement level
 * of any element will change by at most +-1.
*/

#include <t8.h>                            /** General t8code header. Always include this. */
#include <mesh_handle/mesh.hxx>            /** General Mesh header. Always needed for mesh_handle code. */
#include <mesh_handle/competence_pack.hxx> /** Competence Pack for basic mesh_handle features. Look into tutorials/mesh_handle/t8_mesh_competences for more information. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic Cmesh to mesh_handle conversions. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <mesh_handle/concepts.hxx>
#include <t8_types/t8_vec.hxx>          /** t8 vector dataclass. */
#include "t8_mesh_tutorials_common.hxx" /** Default adaption function. */
#include <memory>
#include <span>

/** Build our adapted mesh by transferring the adaption parameters and adapting once with our \ref adapt_callback function.
 * \tparam TMeshClass    The mesh handle class.
 * \param sc_MPI_Comm    The MPI Communicator.
 * \param level          The initial uniform refinement level.
 * \returns Unique pointer to the adapted mesh.
 */
template <t8_mesh_handle::T8MeshType TMeshClass>
std::unique_ptr<TMeshClass>
build_mesh (sc_MPI_Comm comm, int level)
{
  /* Generate a hybrid hypercube, made out of cubes, prisms etc. */
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<TMeshClass> (level, comm);
  /* Defining the adaption parameters. */
  struct adapt_data adapt_params = { { 0.5, 0.5, 1.0 }, 0.2, 0.4 };
  mesh->set_balance ();
  mesh->set_partition ();
  /* Adapting once with our adapt_callback function. */
  mesh->set_adapt (
    TMeshClass::template mesh_adapt_callback_wrapper<adapt_data> (default_adapt_callback<TMeshClass>, adapt_params));
  mesh->set_ghost ();
  mesh->commit ();
  return mesh;
}

/** Entry point of the program. */
int
main (int argc, char **argv)
{
  /*Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /*Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* Print a starting message. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Hello, this is the mesh adaptation example of t8code using the mesh handle.\n");
  t8_global_productionf (" [tutorial] In this example we will adapt a mesh in a spherical shape around a given point "
                         "and write the adapted mesh to a vtu file.\n");
  t8_global_productionf (" [tutorial] \n");

  using mesh_type = t8_mesh_handle::mesh<>;

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Creating an adapted mesh.\n");
  t8_global_productionf (" [tutorial] \n");
  /* The initial uniform refinement level. */
  int uniform_level = 3;
  /* Building the mesh. */
  { /** Scope to ensure mesh is deleted properly. */
    auto mesh = build_mesh<mesh_type> (comm, uniform_level);
    /* Write the mesh to a vtu file. */
    t8_global_productionf (" [tutorial] \n");
    t8_global_productionf (" [tutorial] Writing adapted mesh to vtu file: step3_adapted_mesh.vtu\n");
    t8_global_productionf (" [tutorial] \n");
    t8_mesh_handle::write_mesh_to_vtk (*mesh, "step3_adapted_mesh.vtu");
  }
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
