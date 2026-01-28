/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

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

/**
 * \file t8_gtest_adapt.cxx
 * Tests for the adapt routines of mesh handles. 
 * This tests uses the callback and user data of tutorial step 3 as example.
 * The adaptation criterion is to look at the midpoint coordinates of the current element and if
 * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. 
 */
#include <gtest/gtest.h>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/adapt.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

/** Dummy user data taken from tutorial for test purposes. */
struct dummy_user_data
{
  t8_3D_point midpoint;             /**< The midpoint of our sphere. */
  double refine_if_inside_radius;   /**< If an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /**< If an element's center is larger this value, we coarsen its family. */
};

/** Callback function for the mesh handle to decide for refining or coarsening of (a family of) elements.
 * The function header fits the definition of \ref TMesh::adapt_callback_type_with_userdata.
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
adapt_callback_test ([[maybe_unused]] const TMeshClass &mesh,
                     const std::vector<typename TMeshClass::element_class> &elements, const dummy_user_data &user_data)
{
  auto element_centroid = elements[0].get_centroid ();
  double dist = t8_dist<t8_3D_point, t8_3D_point> (element_centroid, user_data.midpoint);
  if (dist < user_data.refine_if_inside_radius) {
    return 1;
  }
  // Check if we got a family and if yes, if we should coarsen.
  if ((elements.size () > 1) && (dist > user_data.coarsen_if_outside_radius)) {
    return -1;
  }
  return 0;
}

/** Adapt callback implementation for a forest.
 * This callback defines the same adaptation rules as \ref adapt_callback_test,
 * but it is used for the forest instead of the mesh handle.
 */
int
forest_adapt_callback_example (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                               [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                               [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                               [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const struct dummy_user_data *adapt_data = (const struct dummy_user_data *) t8_forest_get_user_data (forest);
  t8_3D_point centroid;
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid.data ());
  double dist = t8_dist<t8_3D_point, t8_3D_point> (centroid, adapt_data->midpoint);
  if (dist < adapt_data->refine_if_inside_radius) {
    return 1;
  }
  else if (is_family && dist > adapt_data->coarsen_if_outside_radius) {
    return -1;
  }
  return 0;
}

/** Test the adapt routine of a mesh handle. 
 * We compare the result to a classically adapted forest with similar callback.
 */
TEST (t8_gtest_handle_adapt, compare_adapt_with_forest)
{
  // Define forest, a mesh handle and user data.
  const int level = 3;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *init_scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, init_scheme, level, 0, sc_MPI_COMM_WORLD);
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<>, dummy_user_data>;
  mesh_class mesh_handle = mesh_class (forest);
  struct dummy_user_data user_data = {
    t8_3D_point ({ 0.5, 0.5, 1 }), /**< Midpoints of the sphere. */
    0.2,                           /**< Refine if inside this radius. */
    0.4                            /**< Coarsen if outside this radius. */
  };

  // Ref the forest as we want to keep using it after the adapt call to compare results.
  t8_forest_ref (forest);

  // Adapt mesh handle.
  mesh_handle.set_adapt (
    mesh_class::mesh_adapt_callback_wrapper<dummy_user_data> (adapt_callback_test<mesh_class>, user_data), false);
  mesh_handle.commit ();
  // Adapt forest classically.
  forest = t8_forest_new_adapt (forest, forest_adapt_callback_example, 0, 1, &user_data);

  // Compare results.
  EXPECT_TRUE (t8_forest_is_equal (mesh_handle.get_forest (), forest));

  // Clean up.
  t8_forest_unref (&forest);
}
