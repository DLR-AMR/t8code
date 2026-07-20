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

/** \file t8_mesh_step3_adapt_forest.cxx
 * This tutorial will explain the the use of competences native to the mesh_handle while explaining the already given competences but also going deeper in how to
 * create custom competences. 
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

using namespace t8_mesh_handle;

/**
 * Creating a simple custom competence that computes the squared volume of an element.
*/
template <typename TUnderlying>
struct volume_squared: public t8_crtp_operator<TUnderlying, volume_squared>
{
 public:
  /**
         * Getter function for the squared volume of an element.
        */
  double
  get_squared_volume () const
  {
    double volume = this->underlying ().get_volume ();
    return volume * volume;
  }
};

/**
 * Demonstrates the use of the standard element data competences by computing the total volume of a mesh.
*/
template <typename MeshType>
void
demonstrate_element_data (MeshType& mesh)
{
  /** The most used geometric standard competences. */
  for (auto& elem : mesh) {
    auto centroid = elem.get_centroid (); /**< Get the Centroid of the element. */
    double volume = elem.get_volume ();   /**< Get the Volume of the element. */

    (void) centroid;

    elem.set_element_data (volume); /**< Save the Volume in the data of the element. */
  }

  double total_volume = 0.0;

  /** Read the element data of each element in the mesh. */
  for (const auto& elem : mesh) {
    total_volume += elem.get_element_data (); /**< Sum up all volumes.*/
  }

  std::cout << "Total volume of the mesh: " << total_volume << '\n';
}

/**
 * Demonstrates the use of the cache competences by comparing the freshly computed values to the one saved in the cache.
*/
template <typename ElementType>
void
demonstrate_cache_competences (const ElementType& elem)
{

  std::cout << "Vertex cache initially filled: " << elem.vertex_cache_filled () << '\n';

  auto vertices1 = elem.get_vertex_coordinates (); /**< Compute the Vertex Coordinates for the first time. */

  std::cout << "Vertex cache after first call: " << elem.vertex_cache_filled () << '\n';

  auto vertices2 = elem.get_vertex_coordinates (); /**< Compute the Vertex Coordinates for the second time. */

  std::cout << "Vertex coordinates (first call):\n";
  for (const auto& v : vertices1) {
    std::cout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
  }

  std::cout << "\nVertex coordinates (second call):\n";
  for (const auto& v : vertices2) {
    std::cout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
  }
}

/**
 * Demonstrates the use of the custom competence 'volume_squared' that was defined at the top so that we can compute the squared volume of each element in the mesh.
*/
template <typename MeshType>
void
demonstrate_custom_competence (const MeshType& mesh)
{
  std::size_t count = 0;

  for (const auto& elem : mesh) {
    if (count++ >= 30) { /**< Counter to only show the first 30 entries so the terminal doesn't get cluttered. */
      break;
    }

    std::cout << "Volume: " << elem.get_volume () /**< Compute default Volume of the element*/
              << "  Squared volume: "
              << elem.get_squared_volume () /**< Computing the squared Volume using the custom competence. */
              << '\n';
  }
}

int
main (int argc, char** argv)
{
  /* The initial refinement level of our meshes. */
  const int level = 2;
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

  /* Print a starting message on the root process. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Hello, this is the competence example of t8code using the mesh handle.\n");
  t8_global_productionf (" [tutorial] In this example we will go through the most important competences and caching,"
                         "as well as create custom competences.\n");
  t8_global_productionf (" [tutorial] \n");

  /* Initializing all the competence packs with the functions/competences we want to use. */
  using cache_competences = element_competence_pack<cache_volume,              /**< Volume Cache. */
                                                    cache_centroid,            /**< Centroid Cache. */
                                                    cache_vertex_coordinates>; /**< Vertex Coordinates Cache. */

  using data_competences = element_competence_pack<
    element_data_element_competence>; /**< Element data Competence to store element data on an element. */

  /** Combine both competence packs into one with union_competence_packs_type. */
  using element_competences = union_competence_packs_type<cache_competences, data_competences>;

  using mesh_competences = mesh_competence_pack<element_data_mesh_competence<
    double>::template type>; /**< Element data Competence to store element data on an element. */

  /* Defining our mesh type with the competence packs defined above. */
  using mesh_type = mesh<element_competences, mesh_competences>;

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Creating a default mesh with refinement level 2.\n");
  t8_global_productionf (" [tutorial] \n");

  /* Creating a simple mesh of Hexahedrons and an initial refinement level of 2. */
  auto default_mesh = handle_hypercube_hybrid_uniform_default<mesh_type> (level, comm, false, false);

  default_mesh->commit (); /**< Committing the mesh so we can work on it. */

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (
    " [tutorial] Demonstrating standard element data competences by computing the total volume.\n");
  t8_global_productionf (" [tutorial] \n");

  demonstrate_element_data (*default_mesh); /**< Calling the element data competence function defined above. */

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Demonstrating the cache competences by comparing the freshly computed values to "
                         "the one saved in the cache.\n");
  t8_global_productionf (" [tutorial] \n");

  demonstrate_cache_competences (*default_mesh->cbegin ()); /**< Calling the cache competence function defined above. */

  /* Defining a competence pack with the volume cache competence and our custom defined competence. */
  using custom_element_competences = element_competence_pack<cache_volume,    /**< Volume cache competence. */
                                                             volume_squared>; /**< Our custom competence. */

  /* Defining a custom mesh_type with our competence pack. */
  using custom_mesh = mesh<custom_element_competences>;

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (
    " [tutorial] Creating a custom mesh for the custom competence with initial refinement level of 2.\n");
  t8_global_productionf (" [tutorial] \n");

  /* Creating a custom mesh with the mesh_type including our custom competence pack and the initial refinement level 2. */
  auto custom = handle_hypercube_hybrid_uniform_default<custom_mesh> (level, comm, false, false);

  custom->commit (); /**< Committing the custom mesh. */

  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Demonstrating the custom competence 'Squared Value'.\n");
  t8_global_productionf (" [tutorial] \n");

  demonstrate_custom_competence (*custom); /**< Calling the custom competence function defined above. */

  /* Finalizing. */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
