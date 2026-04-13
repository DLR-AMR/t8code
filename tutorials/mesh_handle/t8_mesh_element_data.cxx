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
 * This is the same as general/t8_step5_element_data.cxx but using the mesh handle interface instead of the forest 
 * interface.
 */

#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <mesh_handle/mesh_io.hxx>
#include <t8_types/t8_vec.hxx>
#include <memory>
#include <span>

/* The data that we want to store for each element.
 * In this example we want to store the element's level and volume. */
struct data_per_element_type
{
  int level;     /**< Level of the element. */
  double volume; /**< Volume of the element. */
};

/** User data type we will pass to the adapt callback. */
struct user_data
{
  t8_3D_vec midpoint;               /**< The midpoint of our sphere. */
  double refine_if_inside_radius;   /**< If an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /**< If an element's center is larger this value, we coarsen its family. */
};

/** The adaptation callback function. This will refine elements inside of a given sphere and coarsen the elements 
 * outside of a given sphere.
 * \tparam TMeshClass    The mesh handle class.
 * \param [in] mesh      The mesh that should be adapted.
 * \param [in] elements  One element or a family of elements to consider for adaptation.
 * \param [in] user_data The user data to be used during the adaptation process.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *         0 else.
 */
template <typename TMeshClass>
int
adapt_callback ([[maybe_unused]] const TMeshClass &mesh, std::span<const typename TMeshClass::element_class> elements,
                const user_data &user_data)
{
  auto element_centroid = elements[0].get_centroid ();
  double dist = t8_dist<t8_3D_vec, t8_3D_vec> (element_centroid, user_data.midpoint);
  if (dist < user_data.refine_if_inside_radius) {
    return 1;
  }
  // Check if we got a family and if yes, if we should coarsen.
  if ((elements.size () > 1) && (dist > user_data.coarsen_if_outside_radius)) {
    return -1;
  }
  return 0;
}

/** Build a mesh with initial uniform refinement level \a level which is adapted according to \ref adapt_callback, 
 * partitioned and balanced afterwards, and ghost elements are set.
 * \tparam TMeshClass    The mesh handle class.
 * \param [in] comm     MPI communicator to use.
 * \param [in] level Initial refinement level.
 * \return Unique pointer to the mesh created.
 */
template <typename TMeshClass>
std::unique_ptr<TMeshClass>
build_mesh (sc_MPI_Comm comm, int level)
{
  auto mesh_handle = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<TMeshClass> (level, comm);
  struct user_data adapt_data = {
    { 0.5, 0.5, 1 }, /* Midpoint of the sphere. */
    0.2,             /* Refine if inside this radius. */
    0.4              /* Coarsen if outside this radius. */
  };
  /* Adapt, partition, balance and create ghost elements. */
  mesh_handle->set_balance ();
  mesh_handle->set_partition ();
  mesh_handle->set_adapt (
    TMeshClass::template mesh_adapt_callback_wrapper<user_data> (adapt_callback<TMeshClass>, adapt_data));
  mesh_handle->set_ghost ();
  mesh_handle->commit ();
  return mesh_handle;
}

/** Set element data to the mesh handle. 
 * \tparam TMeshClass    The mesh handle class.
 * \param [in, out] mesh  The mesh handle.
 */
template <typename TMeshClass>
void
set_element_data_mesh (TMeshClass &mesh)
{
  for (auto &elem : mesh) {
    elem.set_element_data ({ elem.get_level (), elem.get_volume () });
  }
}

/** Exchange element data set in \ref set_element_data_mesh for ghost elements. 
 * \tparam TMeshClass    The mesh handle class.
 * \param [in, out] mesh  The mesh handle.
 */
template <typename TMeshClass>
void
exchange_ghost_data_mesh (TMeshClass &mesh)
{
  mesh.exchange_ghost_data ();
}

/** Write the mesh as vtu and also write the element's volumes in the file.
 * t8code supports writing element based data to vtu as long as its stored
 * as doubles. Each of the data fields to write has to be provided in its own
 * array of length num_local_elements.
 * We support two types: T8_VTK_SCALAR - One double per element
 *                  and  T8_VTK_VECTOR - 3 doubles per element
 * \tparam TMeshClass     The mesh handle class.
 * \param [in] mesh       The mesh handle.
 * \param [in] fileprefix The prefix of the files where the vtk will be stored.
 *             The master file is then fileprefix.pvtu and the process with rank r writes in the file fileprefix_r.vtu
 */
template <typename TMeshClass>
static void
output_data_to_vtu (const TMeshClass &mesh, const char *prefix)
{
  t8_locidx_t num_elements = mesh.get_num_local_elements ();
  /* We need to allocate a new array to store the volumes on their own.
   * This array has one entry per local element. */
  double *element_volumes = T8_ALLOC (double, num_elements);
  /* The number of user defined data fields to write. */
  int num_data = 1;
  /* For each user defined data field we need one t8_vtk_data_field_t variable */
  t8_vtk_data_field_t vtk_data;
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR. */
  vtk_data.type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data.description, "Element volume");
  vtk_data.data = element_volumes;
  /* Copy the element's volumes from our data array to the output array. */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_volumes[ielem] = mesh[ielem].get_element_data ().volume;
  }
  /* To write user defined data, we need the extended output function write_mesh_to_vtk_ext. 
   * Despite writing user data, it also offers more control over which properties to write. */
  t8_mesh_handle::write_mesh_to_vtk_ext (mesh, prefix, num_data, &vtk_data);
  T8_FREE (element_volumes);
}

/** Entry point of the program. */
int
main (int argc, char **argv)
{
  /* The prefix for our output files. */
  const char *prefix_mesh = "mesh_element_data";
  const char *prefix_mesh_with_data = "mesh_element_data_with_volume_data";
  /* The initial uniform refinement level. */
  const int level = 3;

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

  /* Print a message on the root process. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Hello, this is the element data example of t8code using the mesh handle.\n");
  t8_global_productionf (
    " [tutorial] In this example we will store data on our elements and exchange the data of ghost elements.\n");
  t8_global_productionf (" [tutorial] \n");

  /* Setup: Build cmesh and adapt uniformly. */
  t8_global_productionf (" [tutorial] \n");
  t8_global_productionf (" [tutorial] Creating an adapted mesh.\n");
  t8_global_productionf (" [tutorial] \n");
  { /* We put the mesh in its own scope so that it is automatically destroyed at the end of the scope. 
     * This is only necessary because sc_finalize checks if there are leftover references. 
     * This unique pointer would have been destroyed automatically at the end of the programme. */
    using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<>, data_per_element_type>;
    auto mesh = build_mesh<mesh_class> (comm, level);

    t8_mesh_handle::write_mesh_to_vtk (*mesh, prefix_mesh);
    t8_global_productionf (" [tutorial] Wrote mesh to vtu files: %s*\n", prefix_mesh);

    set_element_data_mesh (*mesh);
    t8_global_productionf (" [tutorial] Computed level and volume data for local elements.\n");
    if (mesh->get_num_local_elements () > 0) {
      /* Output the stored data of the first local element (if it exists). */
      t8_global_productionf (" [tutorial] Element 0 has level %i and volume %e.\n",
                             ((*mesh)[0]).get_element_data ().level, ((*mesh)[0]).get_element_data ().volume);
    }

    /* Exchange the data values of the ghost elements. */
    exchange_ghost_data_mesh (*mesh);
    t8_global_productionf (" [tutorial] Exchanged ghost data.\n");
    if (mesh->get_num_ghosts () > 0) {
      /* Output the data of the first ghost element (if it exists). */
      t8_locidx_t first_ghost_index = mesh->get_num_local_elements ();
      t8_global_productionf (" [tutorial] Ghost 0 has level %i and volume %e.\n",
                             ((*mesh)[first_ghost_index]).get_element_data ().level,
                             ((*mesh)[first_ghost_index]).get_element_data ().volume);
    }

    /* Output the volume data to vtu. */
    output_data_to_vtu (*mesh, prefix_mesh_with_data);
    t8_global_productionf (" [tutorial] Wrote mesh and volume data to %s*.\n", prefix_mesh_with_data);

    /* Cleanup. */
  }  // End scope of mesh
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
