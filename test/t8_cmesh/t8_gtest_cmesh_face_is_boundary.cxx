/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.h>

/* Test for one tree */
class cmesh_face_boundary_one_tree: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();

    /* For each eclass create a cmesh consisting only of one tree. */
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);

    /* We now check each face */
    num_faces = t8_eclass_num_faces[(int) eclass];
  }
  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_eclass eclass;
  int num_faces;
};

TEST_P (cmesh_face_boundary_one_tree, check_face_is_boundary_one_tree)
{

  /* We check whether all faces of the tree are a boundary face. */
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed";

  for (int iface = 0; iface < num_faces; ++iface) {
    ASSERT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 0, iface)) << "Face is not detected as a boundary";
    ASSERT_LT (t8_cmesh_get_face_neighbor (cmesh, 0, iface, NULL, NULL), 0)
      << "Face neighbor on boundary face detected";
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_face_boundary, cmesh_face_boundary_one_tree, AllEclasses, print_eclass);

/* For a two tree cmesh, compute a parallel distribution of the trees.
 * If we have more than 1 process, the first half of process (rank < size/2)
 * gets tree 0, the other half gets tree 1. */
static void
t8_test_compute_parallel_bounds (sc_MPI_Comm comm, t8_gloidx_t *first_tree, t8_gloidx_t *last_tree)
{
  int mpirank;
  int mpisize;
  int mpiret;
  int first_tree_shared;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* If only one process, it gets all trees */
  if (mpisize == 1) {
    *first_tree = 0;
    *last_tree = 1;
    return;
  }
  /* First half of processes gets tree 0, other half gets tree 1 */
  if (mpirank < mpisize / 2) {
    /* The first tree is shared if this rank is in the lower half, but not rank 0 */
    first_tree_shared = mpirank > 0;
    /* The first tree is 0 */
    *first_tree = first_tree_shared ? -1 : 0;
    *last_tree = 0;
  }
  else {
    /* The first tree is shared if this process is not mpisize/2 */
    first_tree_shared = mpirank > mpisize / 2;
    /* The first tree is 1 */
    *first_tree = first_tree_shared ? -2 : 1;
    *last_tree = 1;
  }
}

/* Test for two tree cmesh
 * Creates coarse meshes with two trees for each eclass,
 * one for each face of the first tree as the connecting face.
 * This, only the remaining trees should register as boundary trees. */
class cmesh_face_boundary_two_trees: public testing::TestWithParam<std::tuple<t8_eclass, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    do_partition = std::get<1> (GetParam ());
    num_faces = t8_eclass_num_faces[(int) eclass];
  }

  t8_eclass eclass;
  int do_partition;
  int num_faces;
};

TEST_P (cmesh_face_boundary_two_trees, check_face_is_boundary_two_trees)
{

  t8_cmesh_t cmesh;
  t8_gloidx_t first_tree;
  t8_gloidx_t last_tree;

  t8_test_compute_parallel_bounds (sc_MPI_COMM_WORLD, &first_tree, &last_tree);

  for (int iface = 0; iface < num_faces; ++iface) {
    /* For each face of the eclass we construct one cmesh having this face as a connecting face.
     * Once partitioned and once replicated */
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_tree_class (cmesh, 0, eclass);
    t8_cmesh_set_tree_class (cmesh, 1, eclass);
    /* Connect face iface of tree 0 with face iface of tree 1 with orientation 0 */
    t8_cmesh_set_join (cmesh, 0, 1, iface, iface, 0);
    t8_debugf ("Connecting tree 0 to tree 1 via face %i\n", iface);

    if (do_partition) {
      /* Set the cmesh to be partitioned.
       * We do it in such a way that each process has one local and one ghost tree. */
      t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
    }
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
    ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed";
    for (int checkface = 0; checkface < num_faces; ++checkface) {
      if (iface != checkface) {
        /* The face checkface is a boundary face for tree 0 and tree 1 */
        /* Check that tree 0 face is a boundary */
        ASSERT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 0, checkface)) << "Face is not detected as a boundary";
        /* Check that tree 1 face is a boundary */
        ASSERT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 1, checkface)) << "Face is not detected as a boundary";
        /* Check that we do not detect a face neighbor for tree 0 or tree 1 at this face */
        ASSERT_LT (t8_cmesh_get_face_neighbor (cmesh, 0, checkface, NULL, NULL), 0)
          << "Face neighbor on boundary face detected. Tree 0 face " << checkface << ".";
        ASSERT_LT (t8_cmesh_get_face_neighbor (cmesh, 1, checkface, NULL, NULL), 0)
          << "Face neighbor on boundary face detected. Tree 1 face " << checkface << ".";
      }
      else {
        /* checkface == iface 
         * Thus, tree 0 is connected to tree 1 across this face */
        t8_locidx_t face_neighbor;
        int dual_face = -1, orientation = -1;
        /* Check that tree 0 face is not a boundary */
        ASSERT_FALSE (t8_cmesh_tree_face_is_boundary (cmesh, 0, checkface))
          << "Face is wrongly detected as a boundary.";
        /* Compute the face neighbor info */
        t8_debugf ("Checking face neighbor of local tree 0 across face %i.\n", checkface);
        face_neighbor = t8_cmesh_get_face_neighbor (cmesh, 0, iface, &dual_face, &orientation);
        /* Check the face_neighbor info */
        ASSERT_TRUE (face_neighbor) << "Wrong face neighbor computed. Expected 1 got" << face_neighbor << ".";
        ASSERT_EQ (dual_face, checkface) << "Wrong dual face. Expected " << checkface << " got " << dual_face << ".";
        ASSERT_EQ (orientation, 0) << "Wrong orientation. Expected 0 got " << orientation << ".";
        /* Check that tree 1 face is not a boundary */
        ASSERT_FALSE (t8_cmesh_tree_face_is_boundary (cmesh, 1, checkface))
          << "Face is wrongly detected as a boundary.";
        /* Reset the dual face and orientation to catch false positives (when the get_face_neighbor
         * function does not touch dual_face and orientation) */
        dual_face = orientation = -1;
        /* Compute the face neighbor info */
        t8_debugf ("Checking face neighbor of local tree 1 across face %i.\n", checkface);
        face_neighbor = t8_cmesh_get_face_neighbor (cmesh, 1, checkface, &dual_face, &orientation);
        /* Check the face_neighbor info */
        ASSERT_EQ (face_neighbor, 0) << "Wrong face neighbor computed. Expected 0 got " << face_neighbor << ".";
        ASSERT_EQ (dual_face, checkface) << "Wrong dual face. Expected " << checkface << " got " << dual_face << ".";
        ASSERT_EQ (orientation, 0) << "Wrong orientation. Expected 0 got " << orientation << ".";
      }
    }
    t8_cmesh_destroy (&cmesh);
  } /* End iface loop */
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_face_boundary, cmesh_face_boundary_two_trees,
                          testing::Combine (AllEclasses, testing::Values (0, 1)));
