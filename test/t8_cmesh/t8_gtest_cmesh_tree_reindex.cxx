#include <t8_cmesh/t8_cmesh.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>

#include <gtest/gtest.h>

#include <array>
#include <map>
#include <set>

struct t8_test_cmesh_tree_reindex: public testing::Test
{
 protected:
  static void
  add_tet_tree (t8_cmesh_t cmesh, const t8_gloidx_t global_tree_id, const std::array<double, 12> &vertices)
  {
    t8_cmesh_set_tree_class (cmesh, global_tree_id, T8_ECLASS_TET);
    t8_cmesh_set_tree_vertices (cmesh, global_tree_id, vertices.data (), 4);
  }

  void
  SetUp () override
  {
    t8_cmesh_init (&cmesh);

    ASSERT_NE (cmesh, nullptr);

    /*
     * Build an uncommitted tetrahedral cmesh manually.
     *
     * Each tree has 4 vertices, stored as:
     *
     *   x0, y0, z0, x1, y1, z1, ...
     *
     * The centers are distinct and lie inside the unit cube.
     */
    add_tet_tree (cmesh, 0, { 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0 });

    add_tet_tree (cmesh, 1, { 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 });

    add_tet_tree (cmesh, 2, { 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 });

    add_tet_tree (cmesh, 3, { 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 });

    add_tet_tree (cmesh, 4, { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 });

    add_tet_tree (cmesh, 5, { 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 });
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
};

TEST_F (t8_test_cmesh_tree_reindex, sanity_check_uncommitted_tet_stash)
{
  t8_productionf ("Test started\n");

  constexpr t8_gloidx_t num_trees = 6;

  /*
   * Important:
   * The cmesh must still be uncommitted here.
   * The reindexing function is supposed to read from cmesh->stash.
   */
  EXPECT_FALSE (t8_cmesh_is_committed (cmesh));

  std::map<t8_gloidx_t, t8_gloidx_t> reindex = t8_cmesh_reindex_tree (cmesh, sc_MPI_COMM_SELF);

  /*
   * Reindexing should not commit the original cmesh.
   */
  EXPECT_FALSE (t8_cmesh_is_committed (cmesh));

  /*
   * Every original global tree id should receive exactly one new SFC index.
   */
  EXPECT_EQ (reindex.size (), static_cast<std::size_t> (num_trees));

  std::set<t8_gloidx_t> old_tree_ids;
  std::set<t8_gloidx_t> new_sfc_indices;

  for (const auto &entry : reindex) {
    const t8_gloidx_t old_tree_id = entry.first;
    const t8_gloidx_t new_sfc_index = entry.second;

    EXPECT_GE (old_tree_id, 0);
    EXPECT_LT (old_tree_id, num_trees);

    EXPECT_GE (new_sfc_index, 0);
    EXPECT_LT (new_sfc_index, num_trees);

    old_tree_ids.insert (old_tree_id);
    new_sfc_indices.insert (new_sfc_index);

    t8_productionf ("old global tree id %lli -> new global SFC index %lli\n", static_cast<long long> (old_tree_id),
                    static_cast<long long> (new_sfc_index));
  }

  /*
   * The old tree ids should be exactly {0, 1, 2, 3, 4, 5}.
   */
  EXPECT_EQ (old_tree_ids.size (), static_cast<std::size_t> (num_trees));

  /*
   * The new SFC indices should also be exactly {0, 1, 2, 3, 4, 5}.
   */
  EXPECT_EQ (new_sfc_indices.size (), static_cast<std::size_t> (num_trees));

  for (t8_gloidx_t index = 0; index < num_trees; ++index) {
    EXPECT_TRUE (old_tree_ids.count (index) == 1);
    EXPECT_TRUE (new_sfc_indices.count (index) == 1);
  }
}
