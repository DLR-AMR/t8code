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

/**
 * Tests for each tree class if the boundary conditions input and output remain constant.
 */
TEST_P (t8_cmesh_single_tree_bc, test_single_tree_boundary_conditions)
{
  /* Only check if the tree is local to our process. */
  if (t8_cmesh_get_num_local_trees (cmesh)) {
    /* Retrieve bcs and check that input == output. */
    const auto retrieved_boundary_conditions = t8_cmesh_get_boundary_conditions (cmesh, 0);
    for (size_t i_boundary_condition = 0; i_boundary_condition < boundary_conditions.size (); ++i_boundary_condition) {
      EXPECT_EQ (boundary_conditions[i_boundary_condition], retrieved_boundary_conditions[i_boundary_condition]);
    }
  }
}

/**
 * Tests if internal element faces return empty boundary conditions and if extrior element faces return a filled bc.
 */
TEST_P (t8_cmesh_single_tree_bc, test_single_tree_element_boundary_conditions)
{
  t8_cmesh_ref (cmesh);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 2, 1, sc_MPI_COMM_WORLD);

  /* Only check if the tree is local to our process. */
  if (t8_forest_get_num_local_trees (forest)) {
    t8_locidx_t num_local_elements = t8_forest_get_tree_num_leaf_elements (forest, 0);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);

    /* Some variables for t8_forest_leaf_face_neighbors */
    int *dual_faces;
    int num_neighbors = 0;
    t8_locidx_t *element_indices;
    t8_eclass_t neigh_class;

    for (t8_locidx_t ielem = 0; ielem < num_local_elements; ++ielem) {
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
  }
  t8_forest_unref (&forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_boundary_conditions, t8_cmesh_single_tree_bc, AllEclasses);

/**
 * Tests of a hybrid cmesh with multiple trees if the correct boundary conditions are returned.
 * For the cmesh all faces without a neighbor get the bc "boundary". All faces with a registered neighbor get the bc "internal".
 */
TEST (t8_gtest_cmesh_boundary_conditions, test_hybrid_hypercube_boundary_conditions)
{
  /* We test this with a non periodic, as well as a periodic hypercube. The periodic one should have only "internal" bcs. */
  for (int periodic = 0; periodic < 2; ++periodic) {
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, sc_MPI_COMM_WORLD, periodic);

    /* Test the boundary conditions of the trees. All faces with neighbors should have the bc "internal". All other faces are "boundary". */
    const t8_locidx_t num_local_cmesh_trees = t8_cmesh_get_num_local_trees (cmesh);

    /* Iterate over all cmesh trees. */
    for (t8_locidx_t itree = 0; itree < num_local_cmesh_trees; ++itree) {
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
    t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 1, sc_MPI_COMM_WORLD);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);

    /* Some variables for t8_forest_leaf_face_neighbors */
    int *dual_faces;
    int num_neighbors = 0;
    t8_locidx_t *element_indices;
    t8_eclass_t neigh_class;

    /* Iterate over all trees. */
    const t8_locidx_t num_local_forest_trees = t8_forest_get_num_local_trees (forest);
    for (t8_locidx_t itree = 0; itree < num_local_forest_trees; ++itree) {
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);

      const t8_locidx_t num_local_elements = t8_forest_get_tree_num_leaf_elements (forest, itree);
      /* The forest is level 0, there should only be one element. */
      T8_ASSERT (num_local_elements < 2);

      /* Iterate over our one element. */
      for (t8_locidx_t ielem = 0; ielem < num_local_elements; ++ielem) {
        const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
        const size_t num_faces = scheme->element_get_num_faces (tree_class, elem);

        /* Retrieve the boundary conditions. */
        const auto boundary_conditions = t8_forest_get_boundary_conditions (forest, itree, elem);

        for (size_t iface = 0; iface < num_faces; ++iface) {
          /* Since we have a level 0 forest every face should have a boundary condition. */
          ASSERT_TRUE (boundary_conditions[iface].has_value ());

          /* Find out if we have neighbors. */
          t8_forest_leaf_face_neighbors (forest, itree, elem, NULL, iface, &dual_faces, &num_neighbors,
                                         &element_indices, &neigh_class);
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
    }

    t8_forest_unref (&forest);
  }
}

/**
 * Tests the c interface of the boundary condition module by setting and retrieving the boundary condition of a hex tree.
 * More complex tests are performed for the cpp interface. This test is just to test if th conversion routines are working.
 */
TEST (t8_gtest_cmesh_boundary_conditions, test_boundary_condition_c_interface)
{
  /* Create a cmesh with one hex and apply boundary conditions */
  const char *boundary_conditions[6] = { "bc_0", "bc_1", "bc_2", "bc_3", "bc_4", "bc_5" };
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_boundary_conditions (cmesh, 0, boundary_conditions, 6);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Only check if the cmesh tree is local to our process. */
  if (t8_cmesh_get_num_local_trees (cmesh)) {
    /* Some variables for retrieving and checking the boundary conditions. */
    const char *retrieved_boundary_conditions[6];
    const char *retrieved_single_boundary_condition;
    size_t length = 0;

    /* Retrieve boundary conditions via t8_cmesh_get_boundary_conditions and t8_cmesh_get_boundary_condition() and check them. */
    t8_cmesh_get_boundary_conditions (cmesh, 0, retrieved_boundary_conditions, &length);
    for (size_t i_boundary_condition = 0; i_boundary_condition < length; ++i_boundary_condition) {
      /* Check t8_cmesh_get_boundary_conditions */
      EXPECT_STREQ (boundary_conditions[i_boundary_condition], retrieved_boundary_conditions[i_boundary_condition]);

      /* Check t8_cmesh_get_boundary_condition */
      t8_cmesh_get_boundary_condition (cmesh, 0, i_boundary_condition, &retrieved_single_boundary_condition);
      EXPECT_STREQ (boundary_conditions[i_boundary_condition], retrieved_single_boundary_condition);
    }
  }

  /* We now test the interface for a level 1 forest. This way we should get internal and external faces.
     We assume, that all hex elements inside the tree have the same orientation and face numeration. */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 1, 1, sc_MPI_COMM_WORLD);

  /* Only check if the forest tree is local to our process. */
  if (t8_forest_get_num_local_trees (forest)) {
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    t8_locidx_t num_local_elements = t8_forest_get_tree_num_leaf_elements (forest, 0);

    /* Some variables for t8_forest_leaf_face_neighbors */
    int *dual_faces;
    int num_neighbors = 0;
    t8_locidx_t *element_indices;
    t8_eclass_t neigh_class;

    /* Some variables for retrieving and checking the boundary conditions. */
    const char *retrieved_boundary_conditions[6];
    const char *retrieved_single_boundary_condition;
    size_t length = 0;

    /* Iterate over all elements. */
    for (t8_locidx_t ielem = 0; ielem < num_local_elements; ++ielem) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, 0, ielem);
      const size_t num_faces = scheme->element_get_num_faces (T8_ECLASS_HEX, elem);
      /* Fetch boundary conditions via t8_forest_get_boundary_conditions */
      t8_forest_get_boundary_conditions (forest, 0, elem, retrieved_boundary_conditions, &length);

      for (size_t iface = 0; iface < num_faces; ++iface) {
        t8_forest_leaf_face_neighbors (forest, 0, elem, NULL, iface, &dual_faces, &num_neighbors, &element_indices,
                                       &neigh_class);
        T8_FREE (element_indices);
        T8_FREE (dual_faces);

        /* Fetch boundary conditions via t8_forest_get_boundary_condition */
        t8_forest_get_boundary_condition (forest, 0, elem, iface, &retrieved_single_boundary_condition);

        /* The boundary conditions should be nullptr for internal faces */
        if (num_neighbors > 0) {
          EXPECT_EQ (retrieved_boundary_conditions[iface], nullptr);
          EXPECT_EQ (retrieved_single_boundary_condition, nullptr);
        }
        /* For boundary faces, they should match the bcs of the cmesh cell. */
        else {
          ASSERT_TRUE (retrieved_boundary_conditions[iface] != nullptr);
          ASSERT_TRUE (retrieved_single_boundary_condition != nullptr);

          EXPECT_STREQ (retrieved_boundary_conditions[iface], boundary_conditions[iface]);
          EXPECT_STREQ (retrieved_single_boundary_condition, boundary_conditions[iface]);
        }
      }
    }
  }

  t8_forest_unref (&forest);
}

/**
 * Checks if the registered boundary conditions of the \a handler match with some \a testing_conditions.
 * \param [in]  handler             The handler to check.
 * \param [in]  testing_conditions  The conditions which should be registered in \a handler.
 */
static void
check_boundary_conditions (detail::t8_cmesh_boundary_condition_handler &handler,
                           std::vector<std::string> testing_conditions)
{
  auto retrieved_conditions = handler.get_registered_boundary_conditions ();
  ASSERT_EQ (retrieved_conditions.size (), testing_conditions.size ());

  std::ranges::sort (retrieved_conditions);
  std::ranges::sort (testing_conditions);
  for (size_t i_bc = 0; i_bc < retrieved_conditions.size (); ++i_bc) {
    EXPECT_EQ (retrieved_conditions[i_bc], testing_conditions[i_bc]);
  }
}

/**
 * This test registers boundary conditions on rank 0 and broadcasts them to the other ranks.
 * In the end we check, if every rank hast the boundary conditions.
 */
TEST (t8_gtest_cmesh_boundary_conditions, test_boundary_condition_broadcast)
{
  /* Get MPI info. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int rank;
  sc_MPI_Comm_rank (comm, &rank);

  /* We create the same vector with 20 boundary conditions on each rank. */
  std::vector<std::string> testing_boundary_conditions;
  size_t num_boundary_conditions = 20;
  testing_boundary_conditions.reserve (num_boundary_conditions);
  for (size_t i_bc = 0; i_bc < num_boundary_conditions; ++i_bc) {
    testing_boundary_conditions.emplace_back ("testing_boundary_condition_" + std::to_string (i_bc));
  }

  /* We create a handler and register the boundary conditions only on rank 0. */
  detail::t8_cmesh_boundary_condition_handler handler (nullptr);
  if (rank == 0) {
    for (const auto &boundary_condition : testing_boundary_conditions) {
      handler.register_boundary_condition (boundary_condition);
    }
  }

  /* Broadcast the conditions. */
  handler.bcast (0, comm);

  /* Check if every rank has all conditions. */
  check_boundary_conditions (handler, testing_boundary_conditions);
}

/**
 * This test creates boundary conditions for each rank.
 * Then we synchronize all bcs and check that all boundary conditions are available
 * on all ranks.
 */
TEST (t8_gtest_cmesh_boundary_conditions, test_boundary_condition_synchronize)
{
  /* Get MPI info */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int rank, mpisize;
  sc_MPI_Comm_rank (comm, &rank);
  sc_MPI_Comm_size (comm, &mpisize);

  /* Create the boundary conditions.
   * For each rank, we create 6 bcs.
   * The last rank gets no boundary conditions as a special testing case.
   * Each rank gets his own 6 boundary conditions as well as 2 boundary conditions
   * of the next rank. This way each rank (except the first, last and second to last)
   * has 2 shared bcs with rank - 1, 2 unique bcs and 2 shared bcs with rank + 1.
   * Rank 0 and mpisize - 2 have 4 unique bcs and 2 shared ones (if there are more than 2 ranks).
   * Rank mpisize - 1 has no boundary conditions. */
  std::vector<std::string> testing_boundary_conditions;
  const size_t num_boundary_conditions_per_rank = 6;
  /* Since the last rank gets no bcs we use mpisize - 1 and since every rank gets 2 additional bcs we also add 2. */
  size_t num_boundary_conditions = (mpisize - 1) * num_boundary_conditions_per_rank + 2;
  testing_boundary_conditions.reserve (num_boundary_conditions);
  for (size_t i_bc = 0; i_bc < num_boundary_conditions; ++i_bc) {
    testing_boundary_conditions.emplace_back ("testing_boundary_condition_" + std::to_string (i_bc));
  }

  /* Assign all bcs. */
  detail::t8_cmesh_boundary_condition_handler handler (nullptr);
  size_t bc_min = rank * num_boundary_conditions_per_rank;
  size_t bc_max = (rank + 1) * num_boundary_conditions_per_rank + 2;
  /* All ranks except the last one assign bcs. */
  if (rank != mpisize - 1) {
    for (size_t i_bc = bc_min; i_bc < bc_max; ++i_bc) {
      handler.register_boundary_condition (testing_boundary_conditions[i_bc]);
    }
  }

  /* Synchronize the conditions. */
  handler.synchronize (comm);

  /* Check if every rank has all conditions. */
  check_boundary_conditions (handler, testing_boundary_conditions);
}
