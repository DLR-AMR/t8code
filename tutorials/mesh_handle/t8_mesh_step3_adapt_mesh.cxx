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
 * This is the same as general/t8_step3_adapt_forest.cxx but using the mesh handle interface instead of the forest
 * interface.
*/

#include <t8.h>
#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <mesh_handle/mesh_io.hxx>
#include <mesh_handle/concepts.hxx>
#include <t8_types/t8_vec.hxx>
#include <memory>
#include <span>

/* The data that determines the adaptation characteristics of our algorithm.
 * In this example we want to adapt in a spherical shape around a given point. */
struct adapt_data
{
  t8_3D_vec midpoint;    /**< midpoint of our sphere. */
  double refine_radius;  /**< We refine inside this radius of our sphere.*/
  double coarsen_radius; /**< We coarsen outside this radius of our sphere. */
};

/** The adaption callback function. This will refine elements inside of a given sphere and coarsen the elements 
 * outside of a given sphere.
 * \tparam TMeshClass    The mesh handle class.
 * \param [in] mesh      The mesh that should be adapted.
 * \param [in] elements  One element or a family of elements to consider for adaptation.
 * \param [in] adapt_data The user data to be used during the adaptation process.
 * \returns 1 if the first entry in \a elements should be refined, 
 *         -1 if the family of elements should be coarsened,
 *          0 else.
*/
template <t8_mesh_handle::T8MeshType TMeshClass>
int
adapt_callback ([[maybe_unused]] const TMeshClass &mesh, std::span<const typename TMeshClass::element_class> elements,
                const adapt_data &adapt_data)
{
  auto element_centroid = elements[0].get_centroid ();
  double dist = t8_dist<t8_3D_vec, t8_3D_vec> (element_centroid, adapt_data.midpoint);
  if (dist < adapt_data.refine_radius) {
    return 1; /**< Refine. */
  }           /** First check if there is a family, and only if yes coarsen. */
  else if ((elements.size () > 1) && (dist > adapt_data.coarsen_radius)) {
    return -1; /**< Coarsen. */
  }
  return 0; /**< Do Nothing. */
}

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
  /*Generate a hybrid hypercube, made out of cubes, prisms etc. */
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<TMeshClass> (level, comm);
  /*Defining the adaption parameters. */
  struct adapt_data adapt_params = { { 0.5, 0.5, 1.0 }, 0.2, 0.4 };
  mesh->set_balance ();
  mesh->set_partition ();
  /*Adapting once with our adapt_callback function. */
  mesh->set_adapt (
    TMeshClass::template mesh_adapt_callback_wrapper<adapt_data> (adapt_callback<TMeshClass>, adapt_params));
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
  auto mesh = build_mesh<mesh_type> (comm, uniform_level);
  /* Write the mesh to a vtu file. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Writing adapted mesh to vtu file: step3_adapted_mesh.vtu\n");
  t8_global_productionf (" [tutorial] \n");
  t8_mesh_handle::write_mesh_to_vtk (*mesh, "step3_adapted_mesh.vtu");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
