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

/** \file t8_mesh_step7_interpolation.cxx
 * This is step 7 of the t8code mesh handle tutorials.
 * It is based on general/t8_step7_interpolation.cxx but uses the mesh_handle interface instead. 
 *
 * We showcase how to interpolate/project element data during the adaptation of a mesh.
 * We store a single value on each element (the distance of its centroid to a given point). During refinement
 * the value of the parent is copied to all children; during coarsening the values of the children are averaged
 * and assigned to the parent.
 *
 * Unlike the forest structure, the mesh handle takes care of most things, for example:
 *  - The element data is stored inside the mesh via the element data competences instead of in a user-managed sc_array
 *    (See \ref t8_mesh_stepA_competences.cxx or \ref t8_mesh_element_data.cxx for more information).
 *  - We only need to provide an interpolation callback that receives spans over the old and the new elements.
 *    The mesh calls it automatically on commit, so there is no need to
 *    call t8_forest_iterate_replace ourselves or to compute tree offsets for the data indices.
 */

#include <t8.h>                        /** General t8code header. Always include this. */
#include <mesh_handle/mesh_handle.hxx> /** Convenience header for the mesh_handle. Includes the most necessary headers.  */
#include <t8_types/t8_vec.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <cstring>
#include <span>                         /** Required for the interpolation callback function.  */
#include "t8_mesh_tutorials_common.hxx" /** Adaption function definition used for this tutorial. */

/** The data that we want to store for each element.
 * In this tutorial we store the distance of the element's centroid to a given point.
 * You can store more than one value here (see for example tutorial step 5).
*/
struct step7_element_data_type
{
  double value; /**< Storage for our element value. */
};

/** The mesh class we use in this tutorial.
 *  data_element_competences_basic lets us access the data directly via each element
 *  interpolate_data_mesh_competence_pack stores one element_data_type per element on the mesh and adds the
 *  ability to interpolate this data after adaptation via a user-defined callback. */
using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::data_element_competences_basic,
                                        t8_mesh_handle::interpolate_data_mesh_competence_pack<step7_element_data_type>>;

/** Our interpolation callback function. This decides how element data is recomputed after adaptation. 
 * The mesh calls this function during commit for each group of elements of the old mesh and the corresponding
 * elements of the new adapted mesh.
 * If an element is refined, each child gets the value of its parent.
 * If elements are coarsened, the parent gets the average value of the children.
 * If an element is untouched, it keeps the old value. 
 * For more information on interpolation callbacks, look into \ref t8_mesh_handle::mesh::interpolate_callback_type.
 * \tparam TMeshClass The mesh handle class.
 * \param [in]     mesh_old     The old mesh that is adapted from.
 * \param [in,out] mesh_new     The new mesh constructed from \a mesh_old.
 * \param [in]     refine       -1 if the family got coarsened, 0 if the element was not touched, 1 if it got refined.
 * \param [in]     old_elements Span over the outgoing elements of \a mesh_old
 *                              (the whole family on coarsening, a single element otherwise).
 * \param [in,out] new_elements Span over the incoming elements of \a mesh_new to write the interpolated data to
 *                              (the children on refinement, a single element otherwise).
 */
template <t8_mesh_handle::T8MeshType TMeshClass>
void
interpolate_callback ([[maybe_unused]] const TMeshClass &mesh_old, [[maybe_unused]] TMeshClass &mesh_new,
                      const int refine, std::span<const typename TMeshClass::element_class> old_elements,
                      std::span<typename TMeshClass::element_class> new_elements)
{
  /* Copy old value if the element did not change. */
  if (refine == 0) {
    T8_ASSERT (old_elements.size () == 1 && new_elements.size () == 1);
    new_elements[0].set_element_data (old_elements[0].get_element_data ());
  }
  /* Give the old value to each child if the element was refined. */
  else if (refine == 1) {
    T8_ASSERT (old_elements.size () == 1);
    const step7_element_data_type parent_data = old_elements[0].get_element_data ();
    for (auto &child : new_elements) {
      child.set_element_data (parent_data);
    }
  }
  /* Give the average of all child values to the new element if it was coarsened. */
  else if (refine == -1) {
    T8_ASSERT (new_elements.size () == 1);
    double sum = 0;
    for (const auto &child : old_elements) {
      sum += child.get_element_data ().value;
    }
    new_elements[0].set_element_data (step7_element_data_type { sum / old_elements.size () });
  }
}

/** Compute the initial distance from each element's centroid to a given midpoint. 
 * \tparam TMeshClass    The mesh handle class.
 * \param [in, out] mesh     The mesh handle.
 * \param [in]      midpoint The point to which the distance is computed.
 */
template <t8_mesh_handle::T8MeshType TMeshClass>
void
set_initial_element_data (TMeshClass &mesh, const t8_3D_vec &midpoint)
{
  /* Iterate over each element and set the data directly on the element. */
  for (auto &elem : mesh) {
    const step7_element_data_type data { t8_dist (elem.get_centroid (), midpoint) };
    elem.set_element_data (data);
  }
}

/** Helper function to write the mesh and element data to vtu files. 
 * \tparam TMeshClass   The mesh handle class.
 * \param [in] mesh     The mesh that should be written.
 * \param [in] prefix   Prefix to define the file name.
 */
template <t8_mesh_handle::T8MeshType TMeshClass>
static void
write_vtu (TMeshClass &mesh, const char *prefix)
{
  const t8_locidx_t num_elements = mesh.get_num_local_elements ();

  /* We need to allocate a new array to store the distances on their own.
   * This array has one entry per local element. */
  double *element_values = T8_ALLOC (double, num_elements);
  /* Copy the elements' distances from the mesh to the output array. */
  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    element_values[ielem] = mesh[ielem].get_element_data ().value;
  }

  /* The number of user-defined data fields to write. */
  const int num_data = 1;
  t8_vtk_data_field_t vtk_data; /* For each user data field we need one t8_vtk_data_field_t variable. */
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Centroid distance");
  vtk_data.data = element_values;

  /* To write user-defined data, we need the extended output function write_mesh_to_vtk_ext.
   * Besides writing user data, it also offers more control over which properties to write. */
  const bool write_treeid = true;
  t8_mesh_handle::write_mesh_to_vtk_ext (mesh, prefix, num_data, &vtk_data, write_treeid);

  /* Free the allocated array. */
  T8_FREE (element_values);
}

/** Interpolate the element data after adaptation. As a first step, a uniform mesh is created on a hexahedral
 * coarse mesh. We store the distance of each element's centroid to a given point as element data.
 * The mesh is adapted and the data is interpolated accordingly by the mesh during commit.
 * We write the uniform and the adapted mesh to vtu files.
 */
template <t8_mesh_handle::T8MeshType TMeshClass>
void
t8_step7_mesh_interpolation (TMeshClass &mesh, const char *initial_prefix, const char *adapted_prefix)
{

  /* The parameters for the adaptation. */
  const adapt_user_data adapt_data = { { 0.5, 0.5, 1 }, 0.2, 0.4 };

  /* Compute and set the initial element data. */
  set_initial_element_data (mesh, adapt_data.midpoint);

  /* Write the uniform mesh with its data. */
  write_vtu (mesh, initial_prefix);
  t8_global_productionf (" [step7] Wrote uniform mesh with data to %s*.\n", initial_prefix);

  /* Set the adaptation and interpolation flags. */
  mesh.set_adapt (mesh_class::mesh_adapt_callback_wrapper<adapt_user_data> (adapt_callback<mesh_class>, adapt_data));
  mesh.set_interpolate_callback (interpolate_callback<mesh_class>);
  mesh.set_partition ();

  /* Commit the mesh to change it according to all set flags. */
  mesh.commit ();

  /* Write the adapted mesh with the interpolated data. */
  write_vtu (mesh, adapted_prefix);
  t8_global_productionf (" [step7] Wrote adapted mesh with data to %s*.\n", adapted_prefix);
}

/** Entry point of the program. */
int
main (int argc, char **argv)
{
  const int level = 4; /* Initial refinement level. */
  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  const char *prefix_uniform_mesh = "t8_mesh_step7_uniform_mesh"; /* Prefixes for the vtu files. */
  const char *prefix_adapted_mesh = "t8_mesh_step7_adapted_mesh";

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library; this must happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  { /* Scope for the mesh. It is destroyed and cleaned up automatically at the end of the scope. */
    /* Construct a uniformly refined mesh handle on the coarse mesh using the default scheme. */
    auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, comm);

    /* Define the element data, adapt the mesh and interpolate the element data accordingly. */
    t8_step7_mesh_interpolation (mesh, prefix_uniform_mesh, prefix_adapted_mesh);

  } /* End of mesh scope. */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
