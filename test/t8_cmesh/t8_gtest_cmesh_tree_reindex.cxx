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

INSTANTIATE_TEST_SUITE_P (reordering_modes, t8_test_cmesh_tree_reindex, testing::Values (-1, 0, 1),
                          [] (const testing::TestParamInfo<int> &info) {
                            return info.param < 0 ? "disabled" : (info.param == 0 ? "discard_map" : "retain_map");
                          });
