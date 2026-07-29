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
 * \file t8_gtest_interpolate.cxx
 * Tests to check that the data interpolation for the mesh handle works as intended.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include "t8_gtest_common.hxx"

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <mesh_handle/competences/element_data_competences.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_vec.hxx>
#include <vector>
#include <span>

/** Dummy user data for the interpolation. */
struct interpolate_user_data
{
  int level_step; /**< Levels added when refining and subtracted when coarsening (1 in the standard case). 
                       This is to be applied to the level entry of the dummy user data. */
};

/** Interpolation callback for the mesh handle, using user data.
 * The function header fits the definition of \ref TMesh::interpolate_callback_type_with_userdata.
 * Copies the data of untouched elements, averages the parent volume over the children on refinement, and sums the
 * children's volume onto the parent on coarsening. The level changes by \a user_data.level_step per step.
 * \tparam TMeshClass The mesh handle class.
 * \param [in]     mesh_old     The old mesh that is adapted from.
 * \param [in,out] mesh_new     The new mesh constructed from \a mesh_old.
 * \param [in]     refine       -1 if the family got coarsened, 0 if the element was not touched, 1 if it got refined.
 * \param [in]     old_elements Span over the outgoing elements from \a mesh_old.
 * \param [in,out] new_elements Span over the incoming elements to write the interpolated data to from \a mesh_new.
 * \param [in]     user_data    The user data to be used during the interpolation.
 */
template <typename TMeshClass>
void
interpolate_callback ([[maybe_unused]] const TMeshClass& mesh_old, [[maybe_unused]] TMeshClass& mesh_new,
                      const int refine, std::span<const typename TMeshClass::element_class> old_elements,
                      std::span<typename TMeshClass::element_class> new_elements,
                      const interpolate_user_data& user_data)
{
  /* Untouched: copy data. */
  if (refine == 0) {
    new_elements[0].set_element_data (old_elements[0].get_element_data ());
  }
  /* Refined: children share the parent volume equally, level increases by user_data.level_step. */
  else if (refine == 1) {
    const auto& parent_data = old_elements[0].get_element_data ();
    for (auto& child : new_elements) {
      child.set_element_data (
        data_per_element { parent_data.level + user_data.level_step, parent_data.volume / new_elements.size () });
    }
  }
  /* Coarsened: parent volume is the sum of the children, level decreases by user_data.level_step. */
  else if (refine == -1) {
    double tmp_volume = 0;
    for (const auto& child : old_elements) {
      tmp_volume += child.get_element_data ().volume;
    }
    new_elements[0].set_element_data (
      data_per_element { old_elements[0].get_element_data ().level - user_data.level_step, tmp_volume });
  }
}

TEST (t8_gtest_handle_data, test_interpolate_data)
{
  const int level = 2;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::data_element_competences_basic,
                                          t8_mesh_handle::interpolate_data_mesh_competence_pack<data_per_element>>;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD);

  dummy_user_data user_data_adapt {
    t8_3D_vec { 0.5, 0.5, 1 }, /**< Midpoints of the sphere. */
    0.2,                       /**< Refine if inside this radius. */
    0.4                        /**< Coarsen if outside this radius. */
  };
  interpolate_user_data dummy_interpolate_user_data { 1 }; /**< Level changes by one per adaptation step. */

  // Create element data for all local mesh elements and set via mesh competence.
  std::vector<data_per_element> element_data;
  for (const auto& elem : *mesh) {
    element_data.push_back ({ elem.get_level (), elem.get_volume () });
  }
  mesh->set_element_data (std::move (element_data));

  // Adapt the mesh and set all options.
  mesh->set_adapt (
    mesh_class::mesh_adapt_callback_wrapper<dummy_user_data> (adapt_callback_test<mesh_class>, user_data_adapt));
  mesh->set_balance ();
  mesh->set_interpolate_callback (mesh_class::mesh_interpolate_callback_wrapper<interpolate_user_data> (
    interpolate_callback<mesh_class>, dummy_interpolate_user_data));
  mesh->commit ();

  // Test interpolation.
  bool tested_something = false;
  for (auto& elem : *mesh) {
    if (!tested_something && (elem.get_level () != level)) {
      tested_something = true;
    }
    EXPECT_EQ (elem.get_level (), elem.get_element_data ().level);
    if (!(elem.get_element_data ().level > level)) {
      // For refined elements, the volume is averaged and thus not exact. Only test coarsening.
      EXPECT_EQ (elem.get_volume (), elem.get_element_data ().volume);
    }
  }
  EXPECT_TRUE (tested_something);
}
