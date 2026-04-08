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
 * \file t8_gtest_handle_data.cxx
 * Tests to check that the element data functionality of the \ref t8_mesh_handle::mesh works as intended.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include "t8_gtest_common.hxx"

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <mesh_handle/constructor_wrappers.hxx>
#include <mesh_handle/data_handler.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_types/t8_vec.hxx>
#include <vector>

template <typename TMeshClass>
void
interpolate_callback (const TMeshClass& mesh_old, TMeshClass& mesh_new, const int refine, const int num_old,
                      const t8_locidx_t first_old, const int num_new, const t8_locidx_t first_new)
{
  /* Do not adapt or coarsen. */
  if (refine == 0) {
    mesh_new[first_new].set_element_data (mesh_old[first_old].get_element_data ());
  }
  /* The old element is refined. Volume data is averaged for the children. */
  else if (refine == 1) {
    for (t8_locidx_t i = 0; i < num_new; i++) {
      mesh_new[first_new + i].set_element_data (data_per_element {
        mesh_old[first_old].get_element_data ().level + 1, mesh_old[first_old].get_element_data ().volume / num_new });
    }
  }
  /* Old element is coarsened. */
  else if (refine == -1) {
    double tmp_volume = 0;
    for (t8_locidx_t i = 0; i < num_old; i++) {
      tmp_volume += mesh_old[first_old + i].get_element_data ().volume;
    }
    mesh_new[first_new].set_element_data (
      data_per_element { mesh_old[first_old].get_element_data ().level - 1, tmp_volume });
  }
}

TEST (t8_gtest_handle_data, test_interpolate_data)
{
  const int level = 2;
  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::data_element_competences_basic,
                                          t8_mesh_handle::interpolate_data_mesh_competence<data_per_element>>;
  auto mesh = t8_mesh_handle::handle_hypercube_hybrid_uniform_default<mesh_class> (level, sc_MPI_COMM_WORLD, false,
                                                                                   false, false);

  struct dummy_user_data user_data = {
    t8_3D_vec ({ 0.5, 0.5, 1 }), /**< Midpoints of the sphere. */
    0.2,                         /**< Refine if inside this radius. */
    0.4                          /**< Coarsen if outside this radius. */
  };

  // Create element data for all local mesh elements and set via mesh competence.
  std::vector<data_per_element> element_data;
  for (const auto& elem : *mesh) {
    element_data.push_back ({ elem.get_level (), elem.get_volume () });
  }
  mesh->set_element_data (std::move (element_data));

  mesh->set_adapt (
    mesh_class::mesh_adapt_callback_wrapper<dummy_user_data> (adapt_callback_test<mesh_class>, user_data), false);
  mesh->set_interpolate_callback (interpolate_callback<mesh_class>);
  mesh->commit ();
  bool tested_something = false;
  for (auto& elem : *mesh) {
    if (!tested_something && (elem.get_level () != level)) {
      tested_something = true;
    }
    EXPECT_EQ (elem.get_level (), elem.get_element_data ().level);
    if (!(elem.get_element_data ().level > level)) {
      // For refined elements, the volume is averaged and thus not exact.
      EXPECT_EQ (elem.get_volume (), elem.get_element_data ().volume);
    }
  }
  EXPECT_TRUE (tested_something);
}
