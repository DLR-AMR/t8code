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

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_scheme.hxx>

#include <set>
#include <tuple>
#include <vector>

/**
 * Stores the results of a face iteration.
 */
struct t8_face_test_data
{
  t8_eclass_t eclass;
  int root_face;
  bool stop_recursion;

  int callback_count = 0;
  int leaf_count = 0;

  std::set<t8_locidx_t> leaf_indices;
};

/**
 * Collect visited leaves and verify their validity.
 *
 * \return 1 to continue recursion, 0 to stop.
 */
static int
t8_test_face_callback (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int face,
                       const int is_leaf, const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index,
                       void *user_data)
{
  auto *data = static_cast<t8_face_test_data *> (user_data);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  EXPECT_EQ (ltreeid, 0);
  EXPECT_NE (element, nullptr);
  EXPECT_NE (leaf_elements, nullptr);

  ++data->callback_count;

  /* Every visited element face must belong to the root face. */
  EXPECT_TRUE (scheme->element_is_root_boundary (data->eclass, element, face));

  EXPECT_EQ (scheme->element_get_tree_face (data->eclass, element, face), data->root_face);

  if (is_leaf) {
    EXPECT_GE (tree_leaf_index, 0);
    EXPECT_LT (tree_leaf_index, t8_forest_get_tree_num_leaf_elements (forest, ltreeid));

    ++data->leaf_count;

    /* Each leaf must be visited exactly once. */
    EXPECT_TRUE (data->leaf_indices.insert (tree_leaf_index).second);
  }

  return !data->stop_recursion;
}

/**
 * Parameterized fixture for testing face iteration.
 *
 * Parameters:
 * - Element class
 * - Refinement level
 * - Root face
 * - Stop recursion
 */
class t8_gtest_forest_iterate_faces: public ::testing::TestWithParam<std::tuple<t8_eclass_t, int, int, bool>> {
 protected:
  t8_forest_t forest = nullptr;
  const t8_scheme *scheme = nullptr;
  t8_element_t *root = nullptr;

  t8_eclass_t eclass;

  void
  SetUp () override
  {
    const auto [test_eclass, level, face, stop] = GetParam ();

    eclass = test_eclass;

    /* Construct a single-tree coarse mesh. */
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_class (cmesh, eclass, sc_MPI_COMM_SELF);

    /* Construct a uniformly refined forest. */
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_SELF);

    ASSERT_NE (forest, nullptr);
    ASSERT_EQ (t8_forest_get_num_local_trees (forest), 1);

    scheme = t8_forest_get_scheme (forest);

    /* Construct the root element. */
    scheme->element_new (eclass, 1, &root);
    scheme->set_to_root (eclass, root);

    ASSERT_LT (face, scheme->element_get_num_faces (eclass, root));
  }

  void
  TearDown () override
  {
    if (root != nullptr) {
      scheme->element_destroy (eclass, 1, &root);
    }

    if (forest != nullptr) {
      t8_forest_unref (&forest);
    }
  }

  /**
   * Independently collect the tree-local indices of leaves
   * that touch the specified root face.
   */
  std::set<t8_locidx_t>
  get_expected_leaf_indices (const int root_face) const
  {
    std::set<t8_locidx_t> expected;

    const auto *leaves = t8_forest_tree_get_leaf_elements (forest, 0);
    const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest, 0);

    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf) {
      const t8_element_t *leaf = t8_element_array_index_locidx (leaves, ileaf);

      const int num_faces = scheme->element_get_num_faces (eclass, leaf);

      for (int iface = 0; iface < num_faces; ++iface) {
        if (scheme->element_is_root_boundary (eclass, leaf, iface)
            && scheme->element_get_tree_face (eclass, leaf, iface) == root_face) {
          expected.insert (ileaf);
          break;
        }
      }
    }

    return expected;
  }
};

/**
 * Verify face traversal for all parameter combinations.
 */
TEST_P (t8_gtest_forest_iterate_faces, iterate_faces)
{
  const auto [test_eclass, level, face, stop_recursion] = GetParam ();

  const auto *leaves = t8_forest_tree_get_leaf_elements (forest, 0);

  t8_face_test_data data { test_eclass, face, stop_recursion, 0, 0, {} };

  /* Compute the expected face leaves independently. */
  const auto expected_indices = get_expected_leaf_indices (face);

  /* Execute the face iteration. */
  t8_forest_iterate_faces (forest, 0, root, face, leaves, 0, t8_test_face_callback, &data);

  if (stop_recursion) {
    /* Returning zero must stop recursion at the root. */
    EXPECT_EQ (data.callback_count, 1);
    EXPECT_EQ (data.leaf_count, level == 0 ? 1 : 0);
  }
  else {
    /* All and only the expected leaves must be visited. */
    EXPECT_EQ (data.leaf_indices, expected_indices);
    EXPECT_EQ (data.leaf_count, expected_indices.size ());

    if (level > 0) {
      EXPECT_GT (data.callback_count, data.leaf_count);
    }
    else {
      EXPECT_EQ (data.callback_count, 1);
    }
  }
}

/**
 * Generate test configurations for multiple element classes,
 * refinement levels, faces, and callback behaviors.
 */
static std::vector<std::tuple<t8_eclass_t, int, int, bool>>
t8_generate_face_test_cases ()
{
  std::vector<std::tuple<t8_eclass_t, int, int, bool>> cases;

  const std::vector<std::pair<t8_eclass_t, int>> eclasses
    = { { T8_ECLASS_LINE, 2 }, { T8_ECLASS_QUAD, 4 }, { T8_ECLASS_HEX, 6 } };

  for (const auto &[eclass, num_faces] : eclasses) {
    for (int level = 0; level <= 3; ++level) {
      for (int face = 0; face < num_faces; ++face) {
        for (const bool stop_recursion : { false, true }) {
          cases.emplace_back (eclass, level, face, stop_recursion);
        }
      }
    }
  }

  return cases;
}

INSTANTIATE_TEST_SUITE_P (AllFaceConfigurations, t8_gtest_forest_iterate_faces,
                          ::testing::ValuesIn (t8_generate_face_test_cases ()));
