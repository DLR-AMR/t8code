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

/** \file t8_mesh_stepA_competences.cxx
 * This is step A of the t8code mesh handle tutorials.
 * After finishing the core t8code features, we will now go into an important feature which is native to the mesh handle.
 * These so called competences are a way to extend the functionality of the mesh handle and its elements.
 * 
 * The competences are organized in different types, depending the functionality. 
 * Element data competences are used to store data in the mesh elements and work with it in different ways
 * Cache competences are used to store data in the mesh elements to work with it more efficiently, e.g. to avoid recomputing the same data multiple times.
 * The keypoint about competences though is, that you can create your own competence packs with all the competences you want to use and then use this pack to create a mesh handle with all the functionality you need.
 * This can be further expanded by creating your own competences and adding them to your competence pack, making the mesh handle really flexible and individual for each use case. 
 * 
 * In this tutorial, we will go through the most important competences and caching, as well as create custom competences.
*/

#include <t8.h> /** General t8code header. Always include this. */

#include <mesh_handle/mesh.hxx>                 /** General Mesh header, always needed for mesh_handle code. */
#include <mesh_handle/competence_pack.hxx>      /** Competence Pack for basic mesh_handle features. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic Cmesh to mesh_handle conversions. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <mesh_handle/concepts.hxx> /** Include this to use c++ concepts related to the mesh handle. This can be used to constraint the template parameters to only allow mesh handle classes. */
#include <t8_types/t8_vec.hxx>      /** t8code vector dataclass. */

using namespace t8_mesh_handle; /** Using the namespace to avoid the t8_mesh_handle:: prefix everywhere and shorten the code. */

/**
 * Creating a simple custom competence that computes the squared volume of an element.
 * 
 * All custom competences follow the same CRTP inheritance pattern:
 * They are templated on the underlying element type and inherit from 
 * t8_crtp_operator<TUnderlying, Competence>. This gives the competence access to the functionality
 * of the underlying element with using this->underlying(), allowing it to extend the element with additional methods.
 * The use of t8_crtp_operator also avoids diamond-shaped inheritance when multiple competences are combined into one pack. 
 * 
 * \tparam TUnderlying The underlying element type that we want to extend with this competence.
*/
template <typename TUnderlying>
struct volume_squared_custom_competence: public t8_crtp_operator<TUnderlying, volume_squared_custom_competence>
{
 public:
  /**
  * Returns the squared volume of the underlying element. 
  */
  double
  get_squared_volume () const
  {
    double volume = this->underlying ().get_volume ();
    return volume * volume;
  }
};

/**
 * Example element data type that stores the volume of an element.
*/
struct element_data_volume
{
  double volume; /**< Volume of the element. */
};

/**
 * Demonstrates the use of the standard element data competences by computing the total volume of a mesh.
 * 
 * \param [in] mesh The mesh to compute the total volume of.
 * \param [in] comm The MPI communicator to use for the reduction of the total volume.
*/
template <typename MeshType>
void
demonstrate_element_data (MeshType& mesh, sc_MPI_Comm comm)
{
  /** The most used geometric standard competences. */
  for (auto& elem : mesh) {
    element_data_volume data { elem.get_volume () }; /**< Get the volume of the element. */
    elem.set_element_data (data);                    /**< Save the volume in the data of the element. */
  }

  double local_volume = 0.0;

  /** Read the element data of each element in the mesh. */
  for (const auto& elem : mesh) {
    local_volume += elem.get_element_data ().volume; /**< Sum up all volumes.*/
  }

  double global_volume = 0.0;

  sc_MPI_Reduce (&local_volume, &global_volume, 1, sc_MPI_DOUBLE, sc_MPI_SUM, 0,
                 comm); /**< Reduce the local volumes to the root process. */

  t8_global_productionf (" [t8 Step A Mesh handle] Total volume of the mesh: %f\n", global_volume);
}

/**
 * Demonstrates the use of the cache competences by comparing the freshly computed values to the one saved in the cache.
 * 
 * \param [in] elem The element to demonstrate the cache competences on.
*/
template <typename ElementType>
void
demonstrate_cache_competences (const ElementType& elem)
{

  t8_global_productionf ("Vertex cache initially filled: %d\n", elem.vertex_cache_filled ());

  auto vertices1 = elem.get_vertex_coordinates (); /**< Compute the Vertex Coordinates for the first time. */

  t8_global_productionf ("Vertex coordinates (first call):\n");
  for (const auto& v : vertices1) {
    t8_global_productionf ("(%f, %f, %f)\n", v[0], v[1], v[2]);
  }

  t8_global_productionf ("Vertex cache filled after first call: %d\n", elem.vertex_cache_filled ());

  auto vertices2 = elem.get_vertex_coordinates (); /**< Compute the Vertex Coordinates for the second time. */

  if (vertices1 == vertices2) {
    t8_global_productionf ("Vertex coordinates are the same for both calls.\n");
  }
}

/**
 * Demonstrates the use of the custom competence 'volume_squared' that was defined at the top so that we can compute the squared volume of each element in the mesh.
 * Only the first and last local elements are printed to avoid excessive output when running with multiple MPI processes. 
 * 
 * \param [in] mesh The mesh to demonstrate the custom competence on.
*/
template <typename MeshType>
void
demonstrate_custom_competence (const MeshType& mesh)
{
  auto first_elem = mesh.cbegin ();  /**< Get the first element of this MPI process. */
  auto last_elem = mesh.cend () - 1; /**< Get the last element of this MPI process. */

  t8_global_productionf (
    "First element: Volume: %f Squared volume: %f\n",
    first_elem->get_volume (),          /**< Compute default Volume of the element*/
    first_elem->get_squared_volume ()); /**< Computing the squared Volume using the custom competence. */

  t8_global_productionf (
    "Last element: Volume: %f Squared volume: %f\n",
    last_elem->get_volume (),          /**< Compute default Volume of the element*/
    last_elem->get_squared_volume ()); /**< Computing the squared Volume using the custom competence. */
}

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

  /* Print a starting message on the root process. */
  t8_global_productionf (" [t8 Step A Mesh handle] \n");
  t8_global_productionf (
    " [t8 Step A Mesh handle] Hello, this is the competence tutorial of t8code using the mesh handle.\n");
  t8_global_productionf (
    " [t8 Step A Mesh handle] In this tutorial we will go through the most important competences and caching,"
    "as well as create custom competences.\n");
  t8_global_productionf (" [t8 Step A Mesh handle] \n");
  { /* Start of mesh scope. */
    /* Initializing all the competence packs with the functions/competences we want to use. */

    using data_competences
      = data_element_competences; /**< Element data Competence to store element data on an element. */

    /** Combine the data competence pack with the predefined 'all_cache_element_competences' (see competence_pack.hxx) pack into one with union_competence_packs_type. */
    using element_competences = union_competence_packs_type<all_cache_element_competences, data_competences>;

    using mesh_competences
      = data_mesh_competences<element_data_volume>; /**< Element data Competence to store element data on an element. */

    /* Defining our mesh type with the competence packs defined above. */
    using mesh_type = mesh<element_competences, mesh_competences>;

    const int level = 2;
    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    t8_global_productionf (" [t8 Step A Mesh handle] Creating a default mesh with refinement level %d.\n", level);
    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    /* Creating a simple mesh of Hexahedrons and an initial refinement level of 2. Our competences get transferred onto the mesh by the mesh type we defined above. */
    auto default_mesh = handle_hypercube_hybrid_uniform_default<mesh_type> (level, comm);

    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    t8_global_productionf (
      " [t8 Step A Mesh handle] Demonstrating standard element data competences by computing the total volume.\n");
    t8_global_productionf (" [t8 Step A Mesh handle] \n");

    demonstrate_element_data (*default_mesh, comm); /**< Calling the element data competence function defined above. */

    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    t8_global_productionf (
      " [t8 Step A Mesh handle] Demonstrating the cache competences by comparing the freshly computed values to "
      "the one saved in the cache.\n");
    t8_global_productionf (" [t8 Step A Mesh handle] \n");

    demonstrate_cache_competences (
      (*default_mesh)[0]); /** Only demonstrating the cache competences for the first element of the mesh*/

    /** 
   * We will now create a second mesh with our custom competence pack that includes the volume competence and our custom defined competence 'volume_squared'.
  */
    /* Defining a competence pack with the volume cache competence and our custom defined competence. */
    using custom_element_competences = element_competence_pack<cache_volume, volume_squared_custom_competence>;

    /* Defining a custom mesh_type with our competence pack. */
    using custom_mesh = mesh<custom_element_competences>;

    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    t8_global_productionf (" [t8 Step A Mesh handle] Creating a custom mesh for the custom competence with initial "
                           "refinement level of %d.\n",
                           level);
    t8_global_productionf (" [t8 Step A Mesh handle] \n");

    /* Creating a custom mesh with the mesh_type including our custom competence pack and the initial refinement level 2. */
    auto custom = handle_hypercube_hybrid_uniform_default<custom_mesh> (level, comm);

    t8_global_productionf (" [t8 Step A Mesh handle] \n");
    t8_global_productionf (" [t8 Step A Mesh handle] Demonstrating the custom competence 'Squared Value'.\n");
    t8_global_productionf (" [t8 Step A Mesh handle] \n");

    demonstrate_custom_competence (*custom);
  } /* End of mesh scope. */
  /* Finalizing. */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
