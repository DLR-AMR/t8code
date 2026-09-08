#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>
#include <t8_vtk/t8_vtk_writer.h>

#include <gtest/gtest.h>

#include <array>
#include <map>
#include <set>

struct t8_test_cmesh_tree_reindex: public testing::Test
{
 protected:
  static constexpr std::size_t num_trees = 6;

  static void
  add_tet_tree (t8_cmesh_t cmesh, const t8_gloidx_t global_tree_id, const std::array<double, 12> &vertices)
  {
    t8_cmesh_set_tree_class (cmesh, global_tree_id, T8_ECLASS_TET);
    t8_cmesh_set_tree_vertices (cmesh, global_tree_id, vertices.data (), 4);
  }

  void
  SetUp () override
  {
    /*
     * Same tetrahedral unit-cube decomposition as t8code's
     * t8_cmesh_new_hypercube for T8_ECLASS_TET.
     *
     * Unit cube vertices:
     *
     *   0 = (0, 0, 0)
     *   1 = (1, 0, 0)
     *   2 = (0, 1, 0)
     *   3 = (1, 1, 0)
     *   4 = (0, 0, 1)
     *   5 = (1, 0, 1)
     *   6 = (0, 1, 1)
     *   7 = (1, 1, 1)
     */
    original_vertices = { /*
       * Tree 0: vertices 0, 1, 5, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0 },

                          /*
       * Tree 1: vertices 0, 3, 1, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0 },

                          /*
       * Tree 2: vertices 0, 2, 3, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0 },

                          /*
       * Tree 3: vertices 0, 6, 2, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0 },

                          /*
       * Tree 4: vertices 0, 4, 6, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

                          /*
       * Tree 5: vertices 0, 5, 4, 7
       */
                          std::array<double, 12> { 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 }
    };

    t8_cmesh_init (&cmesh);
    ASSERT_NE (cmesh, nullptr);

    t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

    for (t8_gloidx_t tree_id = 0; tree_id < static_cast<t8_gloidx_t> (num_trees); ++tree_id) {
      add_tet_tree (cmesh, tree_id, original_vertices[static_cast<std::size_t> (tree_id)]);
    }

    /*
     * Same internal face joins as t8code's tetrahedral hypercube.
     */
    t8_cmesh_set_join (cmesh, 0, 1, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 1, 2, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 2, 3, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 3, 4, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 4, 5, 2, 1, 0);
    t8_cmesh_set_join (cmesh, 5, 0, 2, 1, 0);
  }

  void
  TearDown () override
  {
    if (cmesh != nullptr) {
      t8_cmesh_unref (&cmesh);
      cmesh = nullptr;
    }
  }

  t8_cmesh_t cmesh = nullptr;
  std::array<std::array<double, 12>, num_trees> original_vertices;
};

TEST_F (t8_test_cmesh_tree_reindex, commit_reindexes_trees_successfully_and_correctly)
{
  t8_productionf ("Test started\n");

  /*
   * Build an identical uncommitted cmesh to compute the expected reindexing map.
   * This cmesh is later committed without reindexing and written to VTK as a
   * reference output.
   */
  t8_cmesh_t expected_cmesh = nullptr;
  t8_cmesh_init (&expected_cmesh);

  ASSERT_NE (expected_cmesh, nullptr);
  ASSERT_FALSE (t8_cmesh_is_committed (expected_cmesh));

  t8_cmesh_register_geometry<t8_geometry_linear> (expected_cmesh);

  for (t8_gloidx_t tree_id = 0; tree_id < static_cast<t8_gloidx_t> (num_trees); ++tree_id) {
    add_tet_tree (expected_cmesh, tree_id, original_vertices[static_cast<std::size_t> (tree_id)]);
  }

  /*
   * Same internal face joins as t8code's tetrahedral hypercube.
   */
  t8_cmesh_set_join (expected_cmesh, 0, 1, 2, 1, 0);
  t8_cmesh_set_join (expected_cmesh, 1, 2, 2, 1, 0);
  t8_cmesh_set_join (expected_cmesh, 2, 3, 2, 1, 0);
  t8_cmesh_set_join (expected_cmesh, 3, 4, 2, 1, 0);
  t8_cmesh_set_join (expected_cmesh, 4, 5, 2, 1, 0);
  t8_cmesh_set_join (expected_cmesh, 5, 0, 2, 1, 0);

  const std::map<t8_gloidx_t, t8_gloidx_t> expected_reindex = t8_cmesh_reindex_tree (expected_cmesh, sc_MPI_COMM_SELF);

  ASSERT_FALSE (t8_cmesh_is_committed (expected_cmesh));

  /*
   * Check that the computed reindexing map is a valid bijection.
   */
  ASSERT_EQ (expected_reindex.size (), num_trees);

  std::set<t8_gloidx_t> old_tree_ids;
  std::set<t8_gloidx_t> new_tree_ids;

  bool reindex_is_identity = true;

  for (const auto &entry : expected_reindex) {
    const t8_gloidx_t old_tree_id = entry.first;
    const t8_gloidx_t new_tree_id = entry.second;

    EXPECT_GE (old_tree_id, 0);
    EXPECT_LT (old_tree_id, static_cast<t8_gloidx_t> (num_trees));

    EXPECT_GE (new_tree_id, 0);
    EXPECT_LT (new_tree_id, static_cast<t8_gloidx_t> (num_trees));

    old_tree_ids.insert (old_tree_id);
    new_tree_ids.insert (new_tree_id);

    if (old_tree_id != new_tree_id) {
      reindex_is_identity = false;
    }

    t8_productionf ("Expected reindex: old global tree id %lli -> new global tree id %lli\n",
                    static_cast<long long> (old_tree_id), static_cast<long long> (new_tree_id));
  }

  EXPECT_EQ (old_tree_ids.size (), num_trees);
  EXPECT_EQ (new_tree_ids.size (), num_trees);

  for (t8_gloidx_t tree_id = 0; tree_id < static_cast<t8_gloidx_t> (num_trees); ++tree_id) {
    EXPECT_EQ (old_tree_ids.count (tree_id), 1);
    EXPECT_EQ (new_tree_ids.count (tree_id), 1);
  }

  t8_productionf ("Reindexing is %s\n", reindex_is_identity ? "identity" : "non-identity");

  /*
   * Commit the reference cmesh without reindexing and write it to VTK.
   */
  expected_cmesh->reindex_trees = 0;

  t8_cmesh_commit (expected_cmesh, sc_MPI_COMM_SELF);

  ASSERT_TRUE (t8_cmesh_is_committed (expected_cmesh));

  t8_cmesh_vtk_write_file (expected_cmesh, "test_cmesh_tree_reindex_original");

  /*
   * Now commit the actual test cmesh with reindexing enabled.
   */
  ASSERT_FALSE (t8_cmesh_is_committed (cmesh));
  ASSERT_NE (cmesh->stash, nullptr);

  cmesh->reindex_trees = 1;

  t8_cmesh_commit (cmesh, sc_MPI_COMM_SELF);

  ASSERT_TRUE (t8_cmesh_is_committed (cmesh));

  t8_cmesh_vtk_write_file (cmesh, "test_cmesh_tree_reindex_reindexed");

  for (const auto &entry : expected_reindex) {
    const t8_gloidx_t old_tree_id = entry.first;
    const t8_gloidx_t new_tree_id = entry.second;

    const t8_locidx_t new_local_tree_id = t8_cmesh_get_local_id (cmesh, new_tree_id);

    ASSERT_GE (new_local_tree_id, 0);

    EXPECT_EQ (t8_cmesh_get_global_id (cmesh, new_local_tree_id), new_tree_id);

    EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, new_local_tree_id), T8_ECLASS_TET);

    double *actual_vertices = t8_cmesh_get_tree_vertices (cmesh, new_local_tree_id);

    ASSERT_NE (actual_vertices, nullptr);

    const std::array<double, 12> &expected_vertices = original_vertices[static_cast<std::size_t> (old_tree_id)];

    for (int icoord = 0; icoord < 12; ++icoord) {
      EXPECT_DOUBLE_EQ (actual_vertices[icoord], expected_vertices[icoord])
        << "Mismatch for old tree id " << old_tree_id << ", new tree id " << new_tree_id << ", coordinate index "
        << icoord;
    }

    t8_productionf ("Verified old global tree id %lli -> new global tree id %lli\n",
                    static_cast<long long> (old_tree_id), static_cast<long long> (new_tree_id));
  }

  const std::array<std::array<int, 5>, 6> original_joins
    = { std::array<int, 5> { 0, 1, 2, 1, 0 }, std::array<int, 5> { 1, 2, 2, 1, 0 },
        std::array<int, 5> { 2, 3, 2, 1, 0 }, std::array<int, 5> { 3, 4, 2, 1, 0 },
        std::array<int, 5> { 4, 5, 2, 1, 0 }, std::array<int, 5> { 5, 0, 2, 1, 0 } };

  for (const auto &join : original_joins) {
    const t8_gloidx_t old_tree_1 = static_cast<t8_gloidx_t> (join[0]);
    const t8_gloidx_t old_tree_2 = static_cast<t8_gloidx_t> (join[1]);

    const int face_1 = join[2];
    const int face_2 = join[3];
    const int expected_orientation = join[4];

    const t8_gloidx_t new_tree_1 = expected_reindex.at (old_tree_1);
    const t8_gloidx_t new_tree_2 = expected_reindex.at (old_tree_2);

    const t8_locidx_t local_tree_1 = t8_cmesh_get_local_id (cmesh, new_tree_1);

    const t8_locidx_t local_tree_2 = t8_cmesh_get_local_id (cmesh, new_tree_2);

    ASSERT_GE (local_tree_1, 0);
    ASSERT_GE (local_tree_2, 0);

    int dual_face = -1;
    int orientation = -1;

    const t8_locidx_t neighbor_of_tree_1
      = t8_cmesh_get_face_neighbor (cmesh, local_tree_1, face_1, &dual_face, &orientation);

    ASSERT_GE (neighbor_of_tree_1, 0);

    EXPECT_EQ (t8_cmesh_get_global_id (cmesh, neighbor_of_tree_1), new_tree_2);
    EXPECT_EQ (dual_face, face_2);
    EXPECT_EQ (orientation, expected_orientation);

    dual_face = -1;
    orientation = -1;

    const t8_locidx_t neighbor_of_tree_2
      = t8_cmesh_get_face_neighbor (cmesh, local_tree_2, face_2, &dual_face, &orientation);

    ASSERT_GE (neighbor_of_tree_2, 0);

    EXPECT_EQ (t8_cmesh_get_global_id (cmesh, neighbor_of_tree_2), new_tree_1);
    EXPECT_EQ (dual_face, face_1);

    t8_productionf ("Verified reindexed join old trees %lli-%lli -> new trees %lli-%lli\n",
                    static_cast<long long> (old_tree_1), static_cast<long long> (old_tree_2),
                    static_cast<long long> (new_tree_1), static_cast<long long> (new_tree_2));
  }

  t8_cmesh_unref (&expected_cmesh);
}
