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
 * \file t8_gtest_user_data.cxx
 * TODO.
 */

#include <gtest/gtest.h>
#include <t8.h>

#include <mesh_handle/mesh.hxx>
#include <mesh_handle/competence_pack.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>

// Dummy user data taken from a tutorial for test purposes.
struct dummy_user_data
{
  t8_3D_point midpoint;             /* The midpoint of our sphere. */
  double refine_if_inside_radius;   /* if an element's center is smaller than this value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this value, we coarsen its family. */
};

/** Tests that the functionality of the handle gives the same results as if worked with the forest directly.
 * Therefore, we compare the results of the mesh handle with the one accessed with the forest directly.
 */
TEST (t8_gtest_user_data, set_and_get_user_data)
{
  // Define forest and mesh handle.
  const int level = 2;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *init_scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, init_scheme, level, 0, sc_MPI_COMM_WORLD);

  using mesh_class = t8_mesh_handle::mesh<t8_mesh_handle::competence_pack<>, dummy_user_data>;
  mesh_class mesh = mesh_class (forest);

  struct dummy_user_data user_data = {
    t8_3D_point ({ 41, 42, 43 }), /* Midpoints of the sphere. */
    0.2,                          /* Refine if inside this radius. */
    0.4                           /* Coarsen if outside this radius. */
  };

  mesh.set_user_data (&user_data);
  auto mesh_user_data = mesh.get_user_data ();
  EXPECT_EQ (mesh_user_data.midpoint, user_data.midpoint);
  EXPECT_EQ (mesh_user_data.refine_if_inside_radius, user_data.refine_if_inside_radius);
  EXPECT_EQ (mesh_user_data.coarsen_if_outside_radius, user_data.coarsen_if_outside_radius);

  // Unref the forest.
  t8_forest_unref (&forest);
}
