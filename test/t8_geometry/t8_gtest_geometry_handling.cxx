/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_gtest_geometry_handling.cxx
 * This file contains tests for the geometry handling of t8code.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.hxx>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#if T8_ENABLE_OCC
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#endif
/* In this file we collect tests for t8code's cmesh geometry module.
 * These tests are
 *  - test_geometry.test_geometry_handler_register: Tests the geometry_handler register and find interface.
 *  - test_geometry.cmesh_two_trees_and_geometries: Register two geometries for two trees and check that we can access them.
 *  - test_geometry.cmesh_geometry_unique: Check that we can access the geometry via the tree id if
 *                                   we only use one geometry and did not specify tree ids for it.
 *                                   In this case t8code should automatically associate this geometry to all trees.
 */

TEST (test_geometry, test_geometry_handler_register)
{
  t8_geometry_handler geom_handler;

  t8_debugf ("Testing geometry handler register and get geometry.\n");
  /* Throw every implemented geometry at the handler and let it search for it. */
  std::vector<t8_geometry *> geometries;

  geometries.push_back (geom_handler.register_geometry<t8_geometry_linear> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_zero> ());
#if T8_ENABLE_OCC
  geometries.push_back (geom_handler.register_geometry<t8_geometry_cad> ());
#endif /* T8_ENABLE_OCC */
  geometries.push_back (geom_handler.register_geometry<t8_geometry_analytic> ("analytic_geom"));
  geometries.push_back (geom_handler.register_geometry<t8_geometry_linear_axis_aligned> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_lagrange> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_quadrangulated_disk> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_triangulated_spherical_surface> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_tessellated_spherical_surface> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_cubed_spherical_shell> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_prismed_spherical_shell> ());
  geometries.push_back (geom_handler.register_geometry<t8_geometry_cubed_sphere> ());

  /* Check that we can find the geometries by name. */
  for (auto geom : geometries) {
    auto found_geom = geom_handler.get_geometry (geom->t8_geom_get_name ());
    ASSERT_TRUE (found_geom != NULL) << "Could not find registered geometry.";
    ASSERT_EQ (found_geom->t8_geom_get_name (), geom->t8_geom_get_name ())
      << "Could not find geometry with name " << geom->t8_geom_get_name ();
    /* The hash should also be equal. */
    ASSERT_EQ (found_geom->t8_geom_get_hash (), geom->t8_geom_get_hash ())
      << "Could not find geometry with hash " << geom->t8_geom_get_hash ();
  }

  /* Check that we can find the geometries by hash. */
  for (auto geom : geometries) {
    auto found_geom = geom_handler.get_geometry (geom->t8_geom_get_hash ());
    ASSERT_TRUE (found_geom != NULL) << "Could not find registered geometry.";
    ASSERT_EQ (found_geom->t8_geom_get_hash (), geom->t8_geom_get_hash ())
      << "Could not find geometry with hash " << geom->t8_geom_get_hash ();
    /* The name should also be equal. */
    ASSERT_EQ (found_geom->t8_geom_get_name (), geom->t8_geom_get_name ())
      << "Could not find geometry with name " << geom->t8_geom_get_name ();
  }

  /* Try to find a different geometry via the name. Must return nullptr. */
  std::string random_name ("random_name34823412414");
  auto found_geom = geom_handler.get_geometry ("random_name34823412414");
  ASSERT_TRUE (found_geom == nullptr) << "Found a geometry that should not exist.";

  /* Try to find a different geometry via the hash. Must return nullptr. */
  const t8_geometry_hash random_hash (std::hash<std::string> {}(random_name));
  found_geom = geom_handler.get_geometry (random_hash);
  ASSERT_TRUE (found_geom == nullptr) << "Found a geometry that should not exist.";
}

TEST (test_geometry, cmesh_two_trees_and_geometries)
{
  t8_cmesh_t cmesh;

  t8_debugf ("Testing cmesh tree geometry set/get.\n");

  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  /* We will assign a linear geometry, so we need vertices. */
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_QUAD], 4);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  /* Register the linear geometry and zero geometry to this cmesh. */
  auto linear_geom = t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  auto zero_geom = t8_cmesh_register_geometry<t8_geometry_zero> (cmesh);
  /* Set the id geometry for the trees. */
  t8_cmesh_set_tree_geometry (cmesh, 0, linear_geom);
  t8_cmesh_set_tree_geometry (cmesh, 1, zero_geom);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id. */
  const t8_geometry *found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  /* Hash should be equal. */
  ASSERT_EQ (found_geom->t8_geom_get_hash (), linear_geom->t8_geom_get_hash ())
    << "Could not find linear tree geometry at tree 0.";
  /* Name should also be equal. */
  ASSERT_EQ (found_geom->t8_geom_get_name (), linear_geom->t8_geom_get_name ())
    << "Could not find linear tree geometry at tree 0.";

  found_geom = t8_cmesh_get_tree_geometry (cmesh, 1);
  /* Hash should be equal. */
  ASSERT_EQ (found_geom->t8_geom_get_hash (), zero_geom->t8_geom_get_hash ())
    << "Could not find zero tree geometry at tree 0.";
  /* Name should also be equal. */
  ASSERT_EQ (found_geom->t8_geom_get_name (), zero_geom->t8_geom_get_name ())
    << "Could not find zero tree geometry at tree 0.";

  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, cmesh_geometry_unique)
{
  t8_cmesh_t cmesh;

  t8_debugf ("Testing cmesh tree geometry get with unique geometry.\n");

  /* Build a simple 1 tree cmesh and set geometry for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  /* Register the linear_geometry to this cmesh. */
  auto provided_geom = t8_cmesh_register_geometry<t8_geometry_zero> (cmesh);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id.
   * This must now work even though we did not register the geometry for 
   * this tree. Since we only have one geometry. */
  auto found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_TRUE (found_geom != nullptr) << "Could not find any geometry.";
  ASSERT_EQ (found_geom->t8_geom_get_hash (), provided_geom->t8_geom_get_hash ())
    << "Could not find cmesh tree geometry.";

  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, cmesh_no_geometry)
{
  /* Build a simple 1 tree cmesh with no geometry. */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  // Check that the geometry type for tree 0 is invalid
  EXPECT_EQ (t8_geometry_get_type (cmesh, 0), T8_GEOMETRY_TYPE_INVALID);
  // Check that the geometry hash corresponds to invalid geometry
  const t8_geometry_hash hash = t8_cmesh_get_tree_geom_hash (cmesh, 0);
  EXPECT_TRUE (t8_geometry_hash_is_null (hash));
  // Check that we get nullptr when querying the geometry
  const t8_geometry_c *should_be_null = t8_cmesh_get_tree_geometry (cmesh, 0);
  EXPECT_EQ (should_be_null, nullptr);
  // Output to calm the nerves of people looking at the logfiles.
  t8_global_productionf ("We expect an error message about not writing the vtk file here.\n");
  // Try to write vtk file and expect failure
  EXPECT_FALSE (t8_cmesh_vtk_write_file (cmesh, "cmesh_vtk_file"));
  t8_cmesh_vtk_write_file_via_API (cmesh, "cmesh_vtk_file_API", sc_MPI_COMM_WORLD);

  t8_cmesh_destroy (&cmesh);
}
