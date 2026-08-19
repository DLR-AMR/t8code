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

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions_c_interface.h>

/** \file In this file we test the global cmesh vertex numbers.
 *
 * We build a test cmesh consisting of two coarse triangles joined together
 * and associate global vertex numbers with the cmesh's vertices.
 * This cmesh has 4 global vertices in total.
 *
 * We then perform three tests
 *
 * 1) check_tree_to_vertex
 * Here we test the tree_to_vertex connectivity.
 * That is, given a tree id, we get a list of the global vertices of that tree
 * (in local vertex order) and check whether this list is correct.
 *
 * 2) check_vertex_to_tree
 * Here we test the vertex_to_tree connectivity.
 * Given a global vertex index, the vertex_to_tree connectivity returns a list
 * of pairs (local tree_id, local_vertex_id) of all the local trees and their local
 * vertices that are connected to the global vertex.
 * We check whether this list is correct.
 *
 * 3) check_global_vertex_number
 * We verify that the number of global vertices is 4.
 * We additionally verify that the process local number of global vertices is 4 as well.
 * This is true, since the cmesh is not partitioned.
 *
 * Additionally, t8_test_cmesh_vertex_conn_partitioned is the start of a test
 * suite with partitioned cmesh that is currently disabled and could be enabled and extended
 * when cmesh vertex connectivity supports partitioned cmeshes.
 * Note that the test itself then has to be set to parallel in the CMake file.
 */

/**
 * Test fixture for the cmesh boundary condition module. It applies boundary conditions
 * to single tree cmeshes in accordance to their face number.
 */
struct t8_cmesh_single_tree_bc: public testing::TestWithParam<t8_eclass>
{
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    boundary_conditions = { "bc_0", "bc_1", "bc_2", "bc_3", "bc_4", "bc_5" };
    boundary_conditions.resize (static_cast<size_t> (t8_eclass_num_faces[eclass]));
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_tree_class (cmesh, 0, eclass);
    t8_cmesh_set_boundary_conditions (cmesh, 0, boundary_conditions);
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  t8_boundary_conditions<std::string> boundary_conditions;
  t8_cmesh_t cmesh;
  t8_eclass eclass;
};

TEST_P (t8_cmesh_single_tree_bc, test_single_tree_boundary_conditions)
{
  const auto retrieved_boundary_conditions = t8_cmesh_get_boundary_conditions (cmesh, 0);
  for (size_t i_boundary_condition = 0; i_boundary_condition < boundary_conditions.size (); ++i_boundary_condition) {
    EXPECT_EQ (boundary_conditions[i_boundary_condition], retrieved_boundary_conditions[i_boundary_condition]);
  }
}

TEST_P (t8_cmesh_single_tree_bc, test_single_tree_element_boundary_conditions)
{
  t8_cmesh_ref (cmesh);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone (), 2, 0, sc_MPI_COMM_WORLD);
  t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);

  /* Some variables for t8_forest_leaf_face_neighbors */
  int *dual_faces;
  int num_neighbors = 0;
  t8_locidx_t *element_indices;
  t8_eclass_t neigh_class;

  for (t8_locidx_t ielem = 0; ielem < num_elements; ++ielem) {
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, 0, ielem);
    const size_t num_faces = scheme->element_get_num_faces (tree_class, elem);
    for (size_t iface = 0; iface < num_faces; ++iface) {
      t8_forest_leaf_face_neighbors (forest, 0, elem, NULL, iface, &dual_faces, &num_neighbors, &element_indices,
                                     &neigh_class);
      T8_FREE (element_indices);
      T8_FREE (dual_faces);
      const auto boundary_condition = t8_forest_get_boundary_condition (forest, 0, elem, iface);
      if (num_neighbors > 0) {
        EXPECT_FALSE (boundary_condition.has_value ());
      }
      else {
        EXPECT_TRUE (boundary_condition.has_value ());
      }
    }
  }

  t8_forest_unref (&forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_conditions, t8_cmesh_single_tree_bc, AllEclasses);

TEST (t8_gtest_cmesh_boundary_conditions, test_hybrid_hypercube_boundary_conditions)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, 0);

  /* Test the boundary conditions of the trees. All faces with neighbors should have the bc "internal". All other faces are "boundary". */
  const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Iterate over all trees. */
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_faces = t8_eclass_num_faces[tree_class];

    /* Iterate over all faces of the tree. */
    for (int iface = 0; iface < num_faces; ++iface) {
      /* Grab the neighbor eclass and the boundary condition. */
      const t8_eclass neigh_class = t8_cmesh_get_tree_face_neighbor_eclass (cmesh, itree, iface);
      const auto boundary_condition = t8_cmesh_get_boundary_condition (cmesh, itree, iface);
      /* If a face is internal, the boundary condition should be "internal". */
      if (neigh_class == T8_ECLASS_INVALID) {
        EXPECT_EQ (boundary_condition, "boundary");
      }
      else {
        EXPECT_EQ (boundary_condition, "internal");
      }
    }
  }

  /* Do the same test with the forest interface. We set the refinement level to 0 so that there are no internal faces inside trees.
     This way every element face will carry a boundary condition. Internal element faces are checked in another test (test_single_tree_element_boundary_conditions). */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_standalone (), 0, 0, sc_MPI_COMM_WORLD);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  /* Some variables for t8_forest_leaf_face_neighbors */
  int *dual_faces;
  int num_neighbors = 0;
  t8_locidx_t *element_indices;
  t8_eclass_t neigh_class;

  /* Iterate over all trees. */
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);

    /* We will not iterate over the elements, because there is only one per tree. */
    const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, 0);
    const size_t num_faces = scheme->element_get_num_faces (tree_class, elem);

    /* Retrieve the boundary conditions. */
    const auto boundary_conditions = t8_forest_get_boundary_conditions (forest, itree, elem);

    for (size_t iface = 0; iface < num_faces; ++iface) {
      /* Since we have a level 0 forest every face should have a boundary condition. */
      ASSERT_TRUE (boundary_conditions[iface].has_value ());

      /* Find out if we have neighbors. */
      t8_forest_leaf_face_neighbors (forest, itree, elem, NULL, iface, &dual_faces, &num_neighbors, &element_indices,
                                     &neigh_class);
      T8_FREE (element_indices);
      T8_FREE (dual_faces);

      /* If we have neighbors, the bc should be "internal". It should be "boundary" otherwise. */
      if (num_neighbors > 0) {
        EXPECT_EQ (boundary_conditions[iface].value (), "internal");
      }
      else {
        EXPECT_EQ (boundary_conditions[iface].value (), "boundary");
      }
    }
  }

  t8_forest_unref (&forest);
}
