/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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
 * \file t8_gtest_handle_data.cxx
 * Tests to check that \ref t8_mesh_handle::mesh user data and element data functionality works as intended.
 */

#include <gtest/gtest.h>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrapper.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>

// --- Test for user data. ---
/** Dummy user data taken from a tutorial for test purposes. */
struct dummy_user_data
{
  t8_3D_point midpoint;             /**< The midpoint of our sphere. */
  double refine_if_inside_radius;   /**< if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /**< if an element's center is larger this value, we coarsen its family. */
};

/** Check that user data can be set and accesses for the handle.
 */
TEST (t8_gtest_handle_data, set_and_get_user_data)
{
  // Define mesh handle.
  const int level = 2;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<>, dummy_user_data>;
  auto mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);

  struct dummy_user_data user_data = {
    t8_3D_point ({ 41, 42, 43 }), /* Midpoints of the sphere. */
    0.2,                          /* Refine if inside this radius. */
    0.4                           /* Coarsen if outside this radius. */
  };

  // Set user data for the mesh handle and check that the getter returns the same data.
  mesh->set_user_data (&user_data);
  auto mesh_user_data = mesh->get_user_data ();
  EXPECT_EQ (mesh_user_data.midpoint, user_data.midpoint);
  EXPECT_EQ (mesh_user_data.refine_if_inside_radius, user_data.refine_if_inside_radius);
  EXPECT_EQ (mesh_user_data.coarsen_if_outside_radius, user_data.coarsen_if_outside_radius);
}

// --- Test for element data. ---
/** Dummy element data taken from a tutorial for test purposes. */
struct data_per_element
{
  int level;
  double volume;
};

/** Check that element data can be set for the handle and 
 * that the getter has exchanged data for the ghosts.
 */
TEST (t8_gtest_handle_data, set_and_get_element_data)
{
  // Define mesh handle.
  const int level = 2;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<>, void, data_per_element>;
  auto mesh
    = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD, true, true, false);
  if ((mesh->get_dimension () > 1) && (mesh->get_num_local_elements () > 1)) {
    // Ensure that we actually test with ghost elements.
    EXPECT_GT (mesh->get_num_ghosts (), 0);
  }
  auto forest = mesh->get_forest ();

  // Create element data for all local mesh elements.
  std::vector<data_per_element> element_data;
  for (const auto &elem : *mesh) {
    element_data.push_back ({ elem.get_level (), elem.get_volume () });
  }
  mesh->set_element_data (element_data);
  // Get element data and check that the data for all elements (including ghosts) is correct.
  auto mesh_element_data = mesh->exchange_ghost_data ();
  for (t8_locidx_t ielem = 0; ielem < mesh->get_num_local_elements () + mesh->get_num_ghosts (); ielem++) {
    EXPECT_EQ (mesh_element_data[ielem].level, level) << "ielem = " << ielem;
    EXPECT_EQ (mesh_element_data[ielem].volume, (*mesh)[ielem].get_volume ()) << "ielem = " << ielem;
  }
  t8_gloidx_t barrier = t8_forest_get_num_global_trees (forest) / 2.0;
  const int newlevel = 42;
  const double newvolume = 42.42;
  for (auto &elem : *mesh) {
    if (t8_forest_global_tree_id (forest, elem.get_local_tree_id ()) < barrier) {
      elem.set_element_data ({ newlevel, newvolume });
    }
  }
  mesh->exchange_ghost_data ();
  // Check for mesh elements with updated data.
  for (auto &elem : *mesh) {
    if (t8_forest_global_tree_id (forest, elem.get_local_tree_id ()) < barrier) {
      EXPECT_EQ (elem.get_element_data ().level, newlevel);
      EXPECT_EQ (elem.get_element_data ().volume, newvolume);
    }
    else {
      EXPECT_EQ (elem.get_element_data ().level, level);
      EXPECT_EQ (elem.get_element_data ().volume, elem.get_volume ());
    }
  }
  // Check fpor ghost elements with updated data.
  for (t8_locidx_t ighost = mesh->get_num_local_elements ();
       ighost < mesh->get_num_local_elements () + mesh->get_num_ghosts (); ighost++) {
    if (t8_forest_ghost_get_global_treeid (
          forest, (*mesh)[ighost].get_local_tree_id () - t8_forest_get_num_local_trees (forest))
        < barrier) {
      EXPECT_EQ ((*mesh)[ighost].get_element_data ().level, newlevel);
      EXPECT_EQ ((*mesh)[ighost].get_element_data ().volume, newvolume);
    }
    else {
      EXPECT_EQ ((*mesh)[ighost].get_element_data ().level, level);
      EXPECT_EQ ((*mesh)[ighost].get_element_data ().volume, (*mesh)[ighost].get_volume ());
    }
  }
}
