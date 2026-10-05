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
 * After finishing the core t8code features, we will now go into an important feature 
 * which is native to the mesh handle.
 * These so-called competences are a way to extend the functionality of the mesh handle and its elements.
 * 
 * The competences are organized in different types, depending on the functionality. 
 * Element data competences let you attach your own user-defined data to each element and work with it, e.g. exchange it between processes for ghost elements.
 * You may have already seen this in action, because the element data competences were already used in step 5 of the mesh handle tutorials.
 * Cache competences store the results of per-element computations, such as the volume, centroid,
 * vertex coordinates or face neighbors, when they are first requested, so they don't have to be
 * recomputed on later calls. 
 * The key point about competences though is, that you can create your own competence packs with all the competences 
 * you want to use and then use this pack to create a mesh handle with all the functionality you need.
 * This can be further expanded by creating your own competences and adding them to your competence pack, 
 * making the mesh handle really flexible and individual for each use case. 
 * 
 * In this tutorial, we will go through the most important competences and caching, as well as create custom competences.
*/

#include <t8.h> /** General t8code header. Always include this. */

#include <mesh_handle/mesh.hxx>            /** General mesh header, always needed for mesh_handle code. */
#include <mesh_handle/competence_pack.hxx> /** Competence pack for basic mesh_handle features. */
#include <mesh_handle/competences/cache_element_competences.hxx> /** All cache related competences. */
#include <mesh_handle/constructor_wrappers.hxx> /** Wrapper for basic cmesh to mesh_handle conversions. */
#include <mesh_handle/mesh_io.hxx>              /** Used to export mesh to vtk files. */
#include <mesh_handle/concepts.hxx> /** Include this to use c++ concepts related to the mesh handle. This can be used to constrain the template parameters to only allow mesh handle classes. */
#include <t8_types/t8_vec.hxx>      /** t8code vector dataclass. */

using namespace t8_mesh_handle; /** Using the namespace to avoid the t8_mesh_handle:: prefix everywhere and shorten the code. */

/**
 * Creating a simple custom competence that computes the squared volume of an element.
 * 
 * All custom competences have to follow the same CRTP inheritance pattern:
 * They are templated on the underlying element type TUnderlying and inherit from 
 * t8_crtp_operator<TUnderlying, Competence>. This gives the competence access to the functionality
 * of the underlying element with using this->underlying(), allowing it to extend the element with additional methods.
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
 * Stays unused in this tutorial, but necessary to demonstrate competence packs.
 */
struct element_data_volume
{
  double volume; /**< Volume of the element. */
};

/**
 * Demonstrates the use of the cache competences by comparing the freshly computed values to the ones saved in the cache.
 * 
 * \tparam TElementType The mesh element type. 
 * \param [in] elem The element to demonstrate the cache competences on.
*/
template <typename TElementType>
void
demonstrate_cache_competences (const TElementType& elem)
{

  t8_global_productionf (" [mesh_stepA] Vertex cache initially filled: %d\n", elem.vertex_cache_filled ());

  auto vertices1 = elem.get_vertex_coordinates (); /** Compute the vertex coordinates for the first time. */

  t8_global_productionf (" [mesh_stepA] Vertex coordinates (first call):\n");
  for (const auto& v : vertices1) {
    t8_global_productionf ("[mesh_stepA] (%f, %f, %f)\n", v[0], v[1], v[2]);
  }

  t8_global_productionf (" [mesh_stepA] Vertex cache filled after first call: %d\n", elem.vertex_cache_filled ());

  auto vertices2 = elem.get_vertex_coordinates (); /** Compute the vertex coordinates for the second time. */

  if (vertices1 == vertices2) {
    t8_global_productionf (" [mesh_stepA] Vertex coordinates are the same for both calls.\n");
  }
}

/**
 * Demonstrates the use of the custom competence 'volume_squared' that was defined at the top so that we can compute the squared volume of each element in the mesh.
 * Only the first and last local elements of the root process are printed to avoid excessive output when running with multiple MPI processes. 
 * 
 * \tparam TMeshClass The mesh class. 
 * \param [in] mesh The mesh to demonstrate the custom competence.
*/
template <T8MeshType TMeshClass>
void
demonstrate_custom_competence (const TMeshClass& mesh)
{
  /* Guard to prevent undefined behaviour if process has no local elements. */
  if (mesh.cbegin () == mesh.cend ()) {
    return;
  }

  auto first_elem = mesh.cbegin ();  /** Get the first element of this MPI process. */
  auto last_elem = mesh.cend () - 1; /** Get the last element of this MPI process. */

  t8_global_productionf (
    " [mesh_stepA] First element: Volume: %.3e Squared volume: %.3e\n",
    first_elem->get_volume (),          /** Compute default volume of the element*/
    first_elem->get_squared_volume ()); /** Computing the squared volume using the custom competence. */

  t8_global_productionf (
    " [mesh_stepA] Last element: Volume: %.3e Squared volume: %.3e\n",
    last_elem->get_volume (),          /** Compute default volume of the element*/
    last_elem->get_squared_volume ()); /** Computing the squared volume using the custom competence. */
}

int
main (int argc, char** argv)
{
  /* The prefix for our output file. */
  const char* prefix_mesh_with_data = "mesh_competences_with_squared_volume";
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
  t8_global_productionf (" [mesh_stepA] \n");
  t8_global_productionf (" [mesh_stepA] Hello, this is the competence tutorial of t8code using the mesh handle.\n");
  t8_global_productionf (" [mesh_stepA] In this tutorial we will cover the most important competences and caching, "
                         "as well as creating custom competences.\n");
  t8_global_productionf (" [mesh_stepA] \n");
  { /* Start of mesh scope. */
    /* Initializing all the competence packs with the functionality/competences we want to use. */

    /** Combine the data competence pack with the predefined 'all_cache_element_competences' (see competence_pack.hxx) pack into one with union_competence_packs_type. 
     *  Because we already showed the element data competences in tutorial step 5, we will not be showing them again here.
     *  Here we only show how to combine it with other competence packs. 
    */
    using element_competences
      = union_competence_packs_type<all_cache_element_competences, data_element_competences_basic>;

    using mesh_competences
      = data_mesh_competences_basic<element_data_volume>; /**< Mesh competence to store element data on an element. */

    /* Defining our mesh type with the competence packs defined above. */
    using mesh_type = mesh<element_competences, mesh_competences>;

    const int level = 2;
    t8_global_productionf (" [mesh_stepA] \n");
    t8_global_productionf (" [mesh_stepA] Creating a default mesh with refinement level %d.\n", level);
    t8_global_productionf (" [mesh_stepA] \n");
    /* Creating a hybrid uniform mesh. Our competences get transferred onto the mesh by the mesh type we defined above. */
    auto default_mesh = handle_hypercube_hybrid_uniform_default<mesh_type> (level, comm);

    /* Guard to protect from undefined behaviour. */
    if (default_mesh->get_num_local_elements () > 0) {
      t8_global_productionf (" [mesh_stepA] \n");
      t8_global_productionf (
        " [mesh_stepA] Demonstrating the cache competences by comparing the freshly computed values to "
        "the ones saved in the cache.\n");
      t8_global_productionf (" [mesh_stepA] \n");

      demonstrate_cache_competences (
        (*default_mesh)[0]); /** Only demonstrating the cache competences for the first element of the mesh*/
    }
    /** 
      * We will now create a second mesh with our custom competence pack that includes the volume competence and our custom defined competence 'volume_squared'.
    */
    /* Defining a competence pack with the volume cache competence and our custom defined competence. */
    using custom_element_competences = element_competence_pack<cache_volume, volume_squared_custom_competence>;

    /* Defining a custom mesh_type with our competence pack. */
    using custom_mesh_class = mesh<custom_element_competences>;

    t8_global_productionf (" [mesh_stepA] \n");
    t8_global_productionf (" [mesh_stepA] Creating a custom mesh for the custom competence with initial "
                           "refinement level of %d.\n",
                           level);
    t8_global_productionf (" [mesh_stepA] \n");

    /* Creating a mesh with the custom_mesh_class including our custom competence pack. */
    auto custom_mesh = handle_hypercube_hybrid_uniform_default<custom_mesh_class> (level, comm);

    t8_global_productionf (" [mesh_stepA] \n");
    t8_global_productionf (" [mesh_stepA] Demonstrating the custom competence 'squared volume'.\n");
    t8_global_productionf (" [mesh_stepA] \n");

    demonstrate_custom_competence (*custom_mesh);

    /* Collect the squared volume of each local element (computed by our custom competence) for the vtk output. */
    double* squared_volumes = T8_ALLOC (double, custom_mesh->get_num_local_elements ());
    int ielem = 0;
    for (const auto& elem : *custom_mesh) {
      squared_volumes[ielem++] = elem.get_squared_volume ();
    }

    /* Wrap the data in a vtk data field, one scalar value per element. */
    t8_vtk_data_field_t vtk_data;
    vtk_data.type = T8_VTK_SCALAR;
    strcpy (vtk_data.description, "Squared volume");
    vtk_data.data = squared_volumes;

    /* Output the custom competence data to vtu. */
    write_mesh_to_vtk_ext (*custom_mesh, prefix_mesh_with_data, 1, &vtk_data);
    t8_global_productionf (" [mesh_stepA] Wrote mesh and squared volume data to %s*.\n", prefix_mesh_with_data);

    /* Free the memory. */
    T8_FREE (squared_volumes);
  } /* End of mesh scope. */
  /* Finalizing. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
