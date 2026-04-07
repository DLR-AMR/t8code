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
 * \file t8_gtest_adapt_partition_balance.cxx
 * Tests for the adapt, partition and balance routines of mesh handles. 
 * For the adapt routine, we use the callback and user data of tutorial step 3 as example.
 * The adaptation criterion is to look at the midpoint coordinates of the current element and if
 * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen. 
 * The test compares the results of the mesh handle to a forest adapted with the same criterion and balanced and partitioned similarly.
 * Therefore, the check is based on the assumption that the forest functionality works as intended and is tested elsewhere.
 */
#include <gtest/gtest.h>
#include <t8.h>
#include "t8_gtest_common.hxx"

#include <mesh_handle/mesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_vec.hxx>

/** Test the adapt partition and balance routine of a mesh handle. 
 * The test compares the results of the mesh handle to a forest adapted with the same criterion and balanced and partitioned similarly.
 * Therefore, the check is based on the assumption that the forest functionality works as intended and is tested elsewhere.
 */
TEST (t8_gtest_handle_adapt, compare_adapt_with_forest)
{
  // Define forest, a mesh handle and user data.
  const int level = 3;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *init_scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, init_scheme, level, 0, sc_MPI_COMM_WORLD);
  using mesh_class = t8_mesh_handle::mesh<>;
  mesh_class mesh_handle = mesh_class (forest);
  struct dummy_user_data user_data = {
    t8_3D_vec ({ 0.5, 0.5, 1 }), /**< Midpoints of the sphere. */
    0.2,                         /**< Refine if inside this radius. */
    0.4                          /**< Coarsen if outside this radius. */
  };

  // Ref the forest as we want to keep using it after the adapt call to compare results.
  t8_forest_ref (forest);

  // Adapt mesh handle.
  mesh_handle.set_adapt (
    mesh_class::mesh_adapt_callback_wrapper<dummy_user_data> (adapt_callback_test<mesh_class>, user_data), false);
  mesh_handle.commit ();
  // Adapt forest classically.
  forest = t8_forest_new_adapt (forest, forest_adapt_callback_example, 0, 0, &user_data);

  // Compare results.
  EXPECT_TRUE (t8_forest_is_equal (mesh_handle.get_forest (), forest));

  // Adapt the mesh handle again and apply partition and balance.
  mesh_handle.set_balance ();
  mesh_handle.set_partition ();
  mesh_handle.set_adapt (
    mesh_class::mesh_adapt_callback_wrapper<dummy_user_data> (adapt_callback_test<mesh_class>, user_data), false);
  mesh_handle.commit ();
  EXPECT_TRUE (mesh_handle.is_balanced ());

  // Compare the results again to an appropriate forest.
  t8_forest_t forest_compare;
  t8_forest_init (&forest_compare);
  t8_forest_set_user_data (forest_compare, &user_data);
  t8_forest_set_adapt (forest_compare, forest, forest_adapt_callback_example, false);
  t8_forest_set_partition (forest_compare, NULL, false);
  t8_forest_set_balance (forest_compare, NULL, false);
  t8_forest_commit (forest_compare);
  EXPECT_TRUE (t8_forest_is_equal (mesh_handle.get_forest (), forest_compare));

  // Clean up.
  t8_forest_unref (&forest_compare);
}
