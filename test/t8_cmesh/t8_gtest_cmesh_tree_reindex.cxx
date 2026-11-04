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

/** \file t8_gtest_cmesh_tree_reindex.cxx
 * Test opt-in tree reordering and the optional old-to-new tree-ID map.
 */

#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <test/t8_gtest_memory_macros.hxx>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <vector>

/* -1 keeps the default (disabled), 0 enables reordering without retaining the map,
 * and 1 enables reordering and retains the map. */
struct t8_test_cmesh_tree_reindex: public testing::TestWithParam<int>
{
 protected:
  void
  SetUp () override
  {
    t8_cmesh_init (&cmesh);
    t8_cmesh_init (&control_cmesh);
    if (GetParam () >= 0) {
      t8_cmesh_enable_tree_reordering (cmesh, GetParam ());
    }
  }

  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
    t8_cmesh_unref (&control_cmesh);
  }

  static void
  build_lines (t8_cmesh_t mesh, const std::vector<std::array<double, 2>> &endpoints)
  {
    for (t8_gloidx_t tree_id = 0; tree_id < static_cast<t8_gloidx_t> (endpoints.size ()); ++tree_id) {
      const double vertices[6] = { endpoints[tree_id][0], 0, 0, endpoints[tree_id][1], 0, 0 };
      t8_cmesh_set_tree_class (mesh, tree_id, T8_ECLASS_LINE);
      t8_cmesh_set_tree_vertices (mesh, tree_id, vertices, 2);
      t8_cmesh_set_attribute (mesh, tree_id, t8_testsuite_get_package_id (), 0, &tree_id, sizeof (tree_id), 0);

      /* Join touching line segments, even when their original tree IDs are out of order. */
      for (t8_gloidx_t neighbor = 0; neighbor < static_cast<t8_gloidx_t> (endpoints.size ()); ++neighbor) {
        if (endpoints[tree_id][1] == endpoints[neighbor][0]) {
          t8_cmesh_set_join (mesh, tree_id, neighbor, 1, 0, 0);
        }
      }
    }
    t8_cmesh_register_geometry<t8_geometry_linear> (mesh);
    t8_cmesh_commit (mesh, sc_MPI_COMM_SELF);
  }

  void
  check_mesh (const std::vector<t8_gloidx_t> &expected_reindex)
  {
    ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
    ASSERT_TRUE (t8_cmesh_is_committed (control_cmesh));
    const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (control_cmesh);
    ASSERT_EQ (expected_reindex.size (), static_cast<size_t> (num_trees));
    ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), num_trees);
    ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), num_trees);
    EXPECT_EQ (t8_cmesh_get_tree_reindex_map (control_cmesh), nullptr);

    const t8_gloidx_t *reindex_map = t8_cmesh_get_tree_reindex_map (cmesh);
    if (GetParam () == 1) {
      ASSERT_NE (reindex_map, nullptr);
    }
    else {
      EXPECT_EQ (reindex_map, nullptr);
    }

    /* Check the permutation before using it to look up trees or neighbors. */
    std::vector<bool> seen (num_trees, false);
    for (const t8_gloidx_t new_tree_id : expected_reindex) {
      ASSERT_GE (new_tree_id, 0);
      ASSERT_LT (new_tree_id, num_trees);
      ASSERT_FALSE (seen[new_tree_id]);
      seen[new_tree_id] = true;
    }

    for (t8_gloidx_t old_tree_id = 0; old_tree_id < num_trees; ++old_tree_id) {
      const t8_gloidx_t new_tree_id = expected_reindex[old_tree_id];
      SCOPED_TRACE (testing::Message () << "old tree " << old_tree_id << ", new tree " << new_tree_id);
      if (reindex_map != nullptr) {
        EXPECT_EQ (reindex_map[old_tree_id], new_tree_id);
      }
      const t8_locidx_t old_local_id = t8_cmesh_get_local_id (control_cmesh, old_tree_id);
      const t8_locidx_t new_local_id = t8_cmesh_get_local_id (cmesh, new_tree_id);
      ASSERT_GE (old_local_id, 0);
      ASSERT_GE (new_local_id, 0);
      EXPECT_EQ (t8_cmesh_get_global_id (cmesh, new_local_id), new_tree_id);
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (control_cmesh, old_local_id);
      ASSERT_EQ (t8_cmesh_get_tree_class (cmesh, new_local_id), eclass);

      const double *old_vertices = t8_cmesh_get_tree_vertices (control_cmesh, old_local_id);
      const double *new_vertices = t8_cmesh_get_tree_vertices (cmesh, new_local_id);
      ASSERT_NE (old_vertices, nullptr);
      ASSERT_NE (new_vertices, nullptr);
      for (int icoord = 0; icoord < 3 * t8_eclass_num_vertices[eclass]; ++icoord) {
        EXPECT_EQ (new_vertices[icoord], old_vertices[icoord]);
      }

      /* Check every directed face, including boundaries and reverse joins. */
      for (int face = 0; face < t8_eclass_num_faces[eclass]; ++face) {
        SCOPED_TRACE (testing::Message () << "face " << face);
        int old_dual_face = -1, old_orientation = -1;
        int new_dual_face = -1, new_orientation = -1;
        const t8_locidx_t old_neighbor
          = t8_cmesh_get_face_neighbor (control_cmesh, old_local_id, face, &old_dual_face, &old_orientation);
        const t8_locidx_t new_neighbor
          = t8_cmesh_get_face_neighbor (cmesh, new_local_id, face, &new_dual_face, &new_orientation);
        if (old_neighbor < 0) {
          EXPECT_LT (new_neighbor, 0);
        }
        else {
          ASSERT_GE (new_neighbor, 0);
          const t8_gloidx_t old_neighbor_global = t8_cmesh_get_global_id (control_cmesh, old_neighbor);
          EXPECT_EQ (t8_cmesh_get_global_id (cmesh, new_neighbor), expected_reindex[old_neighbor_global]);
          EXPECT_EQ (new_dual_face, old_dual_face);
          EXPECT_EQ (new_orientation, old_orientation);
        }
      }
    }

    /* Retrieving the map must remain possible after querying the committed mesh. */
    const t8_gloidx_t *retrieved_map = t8_cmesh_get_tree_reindex_map (cmesh);
    ASSERT_EQ (retrieved_map, reindex_map);
    if (retrieved_map != nullptr) {
      for (t8_gloidx_t old_tree_id = 0; old_tree_id < num_trees; ++old_tree_id) {
        EXPECT_EQ (retrieved_map[old_tree_id], expected_reindex[old_tree_id]);
      }
    }
  }

  void
  check_lines (const std::vector<std::array<double, 2>> &endpoints, std::vector<t8_gloidx_t> expected_reindex)
  {
    build_lines (control_cmesh, endpoints);
    build_lines (cmesh, endpoints);
    if (GetParam () < 0) {
      std::iota (expected_reindex.begin (), expected_reindex.end (), 0);
    }
    ASSERT_NO_FATAL_FAILURE (check_mesh (expected_reindex));

    /* User attributes must stay attached to the original tree. */
    for (t8_gloidx_t old_tree_id = 0; old_tree_id < static_cast<t8_gloidx_t> (endpoints.size ()); ++old_tree_id) {
      const t8_locidx_t new_local_id = t8_cmesh_get_local_id (cmesh, expected_reindex[old_tree_id]);
      const auto *original_id = static_cast<const t8_gloidx_t *> (
        t8_cmesh_get_attribute (cmesh, t8_testsuite_get_package_id (), 0, new_local_id));
      ASSERT_NE (original_id, nullptr);
      EXPECT_EQ (*original_id, old_tree_id);
    }
  }

  /** Compare a partitioned hypercube against the same ordering on a replicated mesh.
   * \param [in] eclass    Element class used for both hypercubes.
   * \param [in] periodic Whether to connect opposite domain faces.
   */
  void
  check_partitioned_hypercube (const t8_eclass_t eclass, const int periodic)
  {
    int mpirank;
    int mpisize;
    SC_CHECK_MPI (sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank));
    SC_CHECK_MPI (sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize));

    /* Exercise overlapping input records by supplying the complete stash on every rank. */
    if (GetParam () >= 0) {
      t8_cmesh_enable_tree_reordering (control_cmesh, 1);
    }
    t8_cmesh_new_hypercube (&control_cmesh, eclass, sc_MPI_COMM_SELF, 0, 0, periodic);
    t8_cmesh_new_hypercube (&cmesh, eclass, sc_MPI_COMM_WORLD, 0, 1, periodic);
    ASSERT_TRUE (t8_cmesh_is_partitioned (cmesh));
    const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (control_cmesh);
    const t8_gloidx_t first_tree = mpirank * num_trees / mpisize;
    const t8_locidx_t num_local = (mpirank + 1) * num_trees / mpisize - first_tree;
    check_partitioned_mesh (first_tree, num_local);
  }

  /** Compare local trees and ghosts with a replicated reference mesh.
   * \param [in] first_tree  Expected first local global tree ID.
   * \param [in] num_local   Expected number of local trees.
   */
  void
  check_partitioned_mesh (const t8_gloidx_t first_tree, const t8_locidx_t num_local)
  {
    int mpirank;
    SC_CHECK_MPI (sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank));
    ASSERT_TRUE (t8_cmesh_is_partitioned (cmesh));
    const t8_gloidx_t num_trees = t8_cmesh_get_num_trees (control_cmesh);
    ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), num_trees);
    ASSERT_EQ (t8_cmesh_get_first_treeid (cmesh), first_tree);
    ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), num_local);

    /* Check the full map on every rank, including ranks without local trees. */
    const t8_gloidx_t *expected_map = t8_cmesh_get_tree_reindex_map (control_cmesh);
    const t8_gloidx_t *reindex_map = t8_cmesh_get_tree_reindex_map (cmesh);
    if (GetParam () == 1) {
      ASSERT_NE (expected_map, nullptr);
      ASSERT_NE (reindex_map, nullptr);
    }
    else {
      EXPECT_EQ (reindex_map, nullptr);
    }

    std::vector<bool> expected_ghosts (num_trees, false);
    for (t8_locidx_t itree = 0; itree < num_local + t8_cmesh_get_num_ghosts (cmesh); ++itree) {
      const t8_gloidx_t global_tree_id = t8_cmesh_get_global_id (cmesh, itree);
      ASSERT_GE (global_tree_id, 0);
      ASSERT_LT (global_tree_id, num_trees);
      SCOPED_TRACE (testing::Message () << "rank " << mpirank << ", tree " << global_tree_id);
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (control_cmesh, global_tree_id);
      if (itree >= num_local) {
        EXPECT_EQ (t8_cmesh_get_ghost_class (cmesh, itree - num_local), eclass);
        EXPECT_TRUE (expected_ghosts[global_tree_id]);
        expected_ghosts[global_tree_id] = false;
      }
      else {
        EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, itree), eclass);
      }
      const double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
      const double *expected_vertices = t8_cmesh_get_tree_vertices (control_cmesh, global_tree_id);
      ASSERT_NE (vertices, nullptr);
      ASSERT_NE (expected_vertices, nullptr);
      for (int icoordinate = 0; icoordinate < 3 * t8_eclass_num_vertices[eclass]; ++icoordinate) {
        EXPECT_EQ (vertices[icoordinate], expected_vertices[icoordinate]);
      }
      for (int iface = 0; iface < t8_eclass_num_faces[eclass]; ++iface) {
        int dual_face = -1;
        int orientation = -1;
        int expected_dual_face = -1;
        int expected_orientation = -1;
        const t8_locidx_t reference_neighbor = t8_cmesh_get_face_neighbor (control_cmesh, global_tree_id, iface,
                                                                           &expected_dual_face, &expected_orientation);
        const t8_locidx_t neighbor = t8_cmesh_get_face_neighbor (cmesh, itree, iface, &dual_face, &orientation);
        EXPECT_EQ (t8_cmesh_tree_face_is_boundary (cmesh, itree, iface),
                   t8_cmesh_tree_face_is_boundary (control_cmesh, global_tree_id, iface));
        if (itree < num_local && reference_neighbor >= 0
            && (reference_neighbor < first_tree || reference_neighbor >= first_tree + num_local)) {
          expected_ghosts[reference_neighbor] = true;
        }
        const t8_locidx_t expected_neighbor
          = reference_neighbor < 0 ? -1 : t8_cmesh_get_local_id (cmesh, reference_neighbor);
        if (expected_neighbor < 0) {
          EXPECT_LT (neighbor, 0);
        }
        else {
          EXPECT_EQ (neighbor, expected_neighbor);
          EXPECT_EQ (dual_face, expected_dual_face);
          EXPECT_EQ (orientation, expected_orientation);
        }
      }
    }
    EXPECT_EQ (std::count (expected_ghosts.begin (), expected_ghosts.end (), true), 0);
    ASSERT_EQ (t8_cmesh_get_tree_reindex_map (cmesh), reindex_map);
    if (reindex_map != nullptr) {
      ASSERT_NE (expected_map, nullptr);
      for (t8_gloidx_t itree = 0; itree < num_trees; ++itree) {
        EXPECT_EQ (reindex_map[itree], expected_map[itree]);
      }
    }
  }

  /** Check distributed line records, including copied and borrowed variable-size attributes.
   * \param [in] concentrate  Whether only the last rank supplies and owns trees.
   */
  void
  check_distributed_lines (const bool concentrate)
  {
    int mpirank;
    int mpisize;
    SC_CHECK_MPI (sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank));
    SC_CHECK_MPI (sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize));
    const std::vector<std::array<double, 2>> endpoints = { { 3, 4 }, { 0, 1 }, { 5, 6 }, { 1, 2 }, { 4, 5 }, { 2, 3 } };
    std::vector<t8_gloidx_t> expected_reindex = { 3, 0, 5, 1, 4, 2 };
    const t8_gloidx_t num_trees = endpoints.size ();
    const t8_gloidx_t first_tree = concentrate ? 0 : mpirank * num_trees / mpisize;
    const t8_locidx_t num_local
      = concentrate ? (mpirank == mpisize - 1 ? num_trees : 0) : (mpirank + 1) * num_trees / mpisize - first_tree;
    if (GetParam () >= 0) {
      t8_cmesh_enable_tree_reordering (control_cmesh, 1);
    }
    else {
      std::iota (expected_reindex.begin (), expected_reindex.end (), 0);
    }
    build_lines (control_cmesh, endpoints);

    /* Keep only the original local trees and their face ghosts in each input stash. */
    std::vector<bool> supplied_trees (num_trees, false);
    for (t8_gloidx_t itree = first_tree; itree < first_tree + num_local; ++itree) {
      supplied_trees[itree] = true;
      for (t8_gloidx_t ineighbor = 0; ineighbor < num_trees; ++ineighbor) {
        if (endpoints[itree][0] == endpoints[ineighbor][1] || endpoints[itree][1] == endpoints[ineighbor][0]) {
          supplied_trees[ineighbor] = true;
        }
      }
    }
    std::vector<std::vector<unsigned char>> attributes (num_trees);
    t8_cmesh_set_dimension (cmesh, 1);
    for (t8_gloidx_t itree = 0; itree < num_trees; ++itree) {
      attributes[itree].resize (3 * itree + 1);
      for (size_t ibyte = 0; ibyte < attributes[itree].size (); ++ibyte) {
        attributes[itree][ibyte] = static_cast<unsigned char> (17 * itree + ibyte);
      }
      if (supplied_trees[itree]) {
        const double vertices[6] = { endpoints[itree][0], 0, 0, endpoints[itree][1], 0, 0 };
        t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_LINE);
        t8_cmesh_set_tree_vertices (cmesh, itree, vertices, 2);
        t8_cmesh_set_attribute (cmesh, itree, t8_testsuite_get_package_id (), 0, &itree, sizeof (itree), 0);
        t8_cmesh_set_attribute (cmesh, itree, t8_testsuite_get_package_id (), 1, attributes[itree].data (),
                                attributes[itree].size (), itree % 2);
      }
      /* Include outward ghost joins whose remote endpoint may be absent from the local stash. */
      for (t8_gloidx_t ineighbor = 0; ineighbor < num_trees; ++ineighbor) {
        if (endpoints[itree][1] == endpoints[ineighbor][0] && (supplied_trees[itree] || supplied_trees[ineighbor])) {
          t8_cmesh_set_join (cmesh, itree, ineighbor, 1, 0, 0);
        }
      }
    }
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, first_tree + num_local - 1);
    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
    ASSERT_NO_FATAL_FAILURE (check_partitioned_mesh (first_tree, num_local));

    const t8_gloidx_t *reindex_map = t8_cmesh_get_tree_reindex_map (cmesh);
    if (reindex_map != nullptr) {
      for (t8_gloidx_t itree = 0; itree < num_trees; ++itree) {
        EXPECT_EQ (reindex_map[itree], expected_reindex[itree]);
      }
    }
    for (t8_locidx_t itree = 0; itree < num_local + t8_cmesh_get_num_ghosts (cmesh); ++itree) {
      const t8_gloidx_t *original_id
        = static_cast<const t8_gloidx_t *> (t8_cmesh_get_attribute (cmesh, t8_testsuite_get_package_id (), 0, itree));
      ASSERT_NE (original_id, nullptr);
      ASSERT_GE (*original_id, 0);
      ASSERT_LT (*original_id, num_trees);
      EXPECT_EQ (t8_cmesh_get_global_id (cmesh, itree), expected_reindex[*original_id]);
      const unsigned char *attribute
        = static_cast<const unsigned char *> (t8_cmesh_get_attribute (cmesh, t8_testsuite_get_package_id (), 1, itree));
      ASSERT_NE (attribute, nullptr);
      for (size_t ibyte = 0; ibyte < attributes[*original_id].size (); ++ibyte) {
        EXPECT_EQ (attribute[ibyte], attributes[*original_id][ibyte]);
      }
    }
    EXPECT_EQ (t8_cmesh_get_tree_reindex_map (cmesh), reindex_map);
  }

  t8_cmesh_t cmesh = nullptr;
  t8_cmesh_t control_cmesh = nullptr;
};

TEST_P (t8_test_cmesh_tree_reindex, scrambled_lines)
{
  /* The 1D SFC orders these segments from left to right. */
  check_lines ({ { 3, 4 }, { 0, 1 }, { 5, 6 }, { 1, 2 }, { 4, 5 }, { 2, 3 } }, { 3, 0, 5, 1, 4, 2 });
}

TEST_P (t8_test_cmesh_tree_reindex, small_bounding_box)
{
  check_lines ({ { 2e-10, 3e-10 }, { 0, 1e-10 }, { 1e-10, 2e-10 } }, { 2, 0, 1 });
}

TEST_P (t8_test_cmesh_tree_reindex, coincident_centers)
{
  /* Trees 0, 2, and 3 have the same center. Keep all of them in their original ID order. */
  check_lines ({ { 2, 4 }, { 0, 2 }, { 2.5, 3.5 }, { 1, 5 }, { 4, 6 } }, { 1, 0, 2, 3, 4 });
}

TEST_P (t8_test_cmesh_tree_reindex, single_tree)
{
  check_lines ({ { 0, 1 } }, { 0 });
}

TEST_P (t8_test_cmesh_tree_reindex, hypercube_preserves_geometry_and_connectivity)
{
  /* Periodic joins also exercise nonzero face orientations. */
  t8_cmesh_new_hypercube (&control_cmesh, T8_ECLASS_TET, sc_MPI_COMM_SELF, 0, 0, 1);
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_TET, sc_MPI_COMM_SELF, 0, 0, 1);
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), 6);

  /* Reconstruct the actual permutation independently of the retained map. */
  using vertex_array_t = std::array<double, 12>;
  std::map<vertex_array_t, t8_gloidx_t> vertices_to_old_id;
  for (t8_locidx_t old_id = 0; old_id < 6; ++old_id) {
    const double *vertices = t8_cmesh_get_tree_vertices (control_cmesh, old_id);
    ASSERT_NE (vertices, nullptr);
    vertex_array_t key;
    std::copy_n (vertices, key.size (), key.begin ());
    ASSERT_TRUE (vertices_to_old_id.emplace (key, old_id).second);
  }
  std::vector<t8_gloidx_t> observed_reindex (6, -1);
  bool changed = false;
  for (t8_locidx_t new_id = 0; new_id < 6; ++new_id) {
    const double *vertices = t8_cmesh_get_tree_vertices (cmesh, new_id);
    ASSERT_NE (vertices, nullptr);
    vertex_array_t key;
    std::copy_n (vertices, key.size (), key.begin ());
    const auto old_tree = vertices_to_old_id.find (key);
    ASSERT_NE (old_tree, vertices_to_old_id.end ());
    observed_reindex[old_tree->second] = new_id;
    changed = changed || old_tree->second != new_id;
  }
  EXPECT_EQ (changed, GetParam () >= 0);
  check_mesh (observed_reindex);
}

TEST_P (t8_test_cmesh_tree_reindex, partitioned_hypercube)
{
  check_partitioned_hypercube (T8_ECLASS_TET, 1);
}

TEST_P (t8_test_cmesh_tree_reindex, partitioned_single_tree)
{
  check_partitioned_hypercube (T8_ECLASS_HEX, 0);
}

TEST_P (t8_test_cmesh_tree_reindex, partitioned_distributed_lines)
{
  check_distributed_lines (false);
}

TEST_P (t8_test_cmesh_tree_reindex, partitioned_lines_empty_root)
{
  check_distributed_lines (true);
}

TEST_P (t8_test_cmesh_tree_reindex, partitioned_rank_local_map_retention)
{
  int mpirank;
  int mpisize;
  SC_CHECK_MPI (sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank));
  SC_CHECK_MPI (sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize));
  const bool retain_map = GetParam () == 1 && mpirank == mpisize - 1;
  if (GetParam () >= 0) {
    t8_cmesh_enable_tree_reordering (cmesh, retain_map ? -1 : 0);
  }
  t8_cmesh_new_hypercube (&cmesh, T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 1, 0);
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));
  EXPECT_EQ (t8_cmesh_get_num_trees (cmesh), 1);
  EXPECT_EQ (t8_cmesh_get_num_local_trees (cmesh), mpirank == mpisize - 1 ? 1 : 0);
  EXPECT_EQ (t8_cmesh_get_num_ghosts (cmesh), 0);
  const t8_gloidx_t *reindex_map = t8_cmesh_get_tree_reindex_map (cmesh);
  if (retain_map) {
    ASSERT_NE (reindex_map, nullptr);
    EXPECT_EQ (reindex_map[0], 0);
  }
  else {
    EXPECT_EQ (reindex_map, nullptr);
  }
  if (mpirank == mpisize - 1) {
    EXPECT_EQ (t8_cmesh_get_global_id (cmesh, 0), 0);
    for (int iface = 0; iface < t8_eclass_num_faces[T8_ECLASS_HEX]; ++iface) {
      EXPECT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 0, iface));
    }
  }
  ASSERT_EQ (t8_cmesh_get_tree_reindex_map (cmesh), reindex_map);
  if (reindex_map != nullptr) {
    EXPECT_EQ (reindex_map[0], 0);
  }
}

INSTANTIATE_TEST_SUITE_P (reordering_modes, t8_test_cmesh_tree_reindex, testing::Values (-1, 0, 1),
                          [] (const testing::TestParamInfo<int> &info) {
                            return info.param < 0 ? "disabled" : (info.param == 0 ? "discard_map" : "retain_map");
                          });
