#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_tree_reindex.hxx>

#include <gtest/gtest.h>

#include <map>
#include <set>

struct t8_test_cmesh_tree_reindex: public testing::Test
{
 protected:
  void
  SetUp () override
  {
    const int do_partition = 0;
    const int periodic = 0;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_SELF, do_partition, periodic, 0);

    ASSERT_NE (cmesh, nullptr);
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

TEST_F (t8_test_cmesh_tree_reindex, sanity_check_tet_hypercube)
{
  t8_productionf ("Test started\n");

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);

  /*
   * A tetrahedral hypercube consists of 6 tetrahedron trees.
   */
  EXPECT_EQ (num_local_trees, 6);

  std::map<t8_locidx_t, t8_locidx_t> reindex = t8_cmesh_reindex_tree (cmesh, sc_MPI_COMM_SELF);

  /*
   * Every original local tree should receive exactly one new SFC index.
   */
  EXPECT_EQ (reindex.size (), static_cast<std::size_t> (num_local_trees));

  std::set<t8_locidx_t> old_tree_ids;
  std::set<t8_locidx_t> new_sfc_indices;

  for (const auto &entry : reindex) {
    const t8_locidx_t old_tree_id = entry.first;
    const t8_locidx_t new_sfc_index = entry.second;

    EXPECT_GE (old_tree_id, 0);
    EXPECT_LT (old_tree_id, num_local_trees);

    EXPECT_GE (new_sfc_index, 0);
    EXPECT_LT (new_sfc_index, num_local_trees);

    old_tree_ids.insert (old_tree_id);
    new_sfc_indices.insert (new_sfc_index);

    t8_productionf ("old local tree id %li -> new local SFC index %li\n", static_cast<long> (old_tree_id),
                    static_cast<long> (new_sfc_index));
  }

  /*
   * The old tree ids should be exactly {0, 1, 2, 3, 4, 5}.
   */
  EXPECT_EQ (old_tree_ids.size (), static_cast<std::size_t> (num_local_trees));

  /*
   * The new SFC indices should also be exactly {0, 1, 2, 3, 4, 5}.
   */
  EXPECT_EQ (new_sfc_indices.size (), static_cast<std::size_t> (num_local_trees));

  for (t8_locidx_t index = 0; index < num_local_trees; ++index) {
    EXPECT_TRUE (old_tree_ids.count (index) == 1);
    EXPECT_TRUE (new_sfc_indices.count (index) == 1);
  }
}
