
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
#include <utility>
#include <vector>

/**
 * Stores the configuration and results of a face iteration.
 *
 * An instance of this structure is passed to the iteration callback
 * through its user_data argument. The callback uses the configuration
 * to validate visited elements and records information about the
 * traversal for subsequent assertions.
 */
struct t8_face_test_data
{
  /** The element class of the tree being traversed. */
  t8_eclass_t eclass;

  /** The root face on which the iteration is performed. */
  int root_face;

  /** Whether the callback should terminate recursion immediately. */
  bool stop_recursion;

  /** The total number of callback invocations, including non-leaf elements. */
  int callback_count = 0;

  /** The number of leaf elements visited during the iteration. */
  int leaf_count = 0;

  /** The unique tree-local indices of all visited leaf elements. */
  std::set<t8_locidx_t> leaf_indices;
};

/**
 * Callback function used to validate the face iteration.
 *
 * This callback is invoked for each element encountered during the
 * traversal of a specified root face. It verifies that the element
 * belongs to the expected tree and that its current face corresponds
 * to the requested root face.
 *
 * For leaf elements, the callback additionally verifies the validity
 * of the tree-local leaf index and records it in a set. This allows
 * the test to detect duplicate visits and compare the visited leaves
 * against an independently computed reference set.
 *
 * The callback can optionally terminate recursion by returning zero.
 * This behavior is used to verify that t8_forest_iterate_faces
 * respects the callback's return value.
 *
 * \param [in]     forest           The forest containing the element.
 * \param [in]     ltreeid          The local index of the tree being traversed.
 * \param [in]     element          The current element visited by the iterator.
 * \param [in]     face             The face of the current element being traversed.
 * \param [in]     is_leaf          Nonzero if the current element is a leaf.
 * \param [in]     leaf_elements    The array of leaf elements associated with
 *                                 the current element.
 * \param [in]     tree_leaf_index  The tree-local index of the first leaf in
 *                                 the current element's leaf array.
 * \param [in,out] user_data        A pointer to the test configuration and
 *                                 traversal results.
 *
 * \return                        Zero if recursion should terminate at the
 *                                current element, nonzero otherwise.
 */
static int
t8_test_face_callback (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int face,
                       const int is_leaf, const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index,
                       void *user_data)
{
  auto *data = static_cast<t8_face_test_data *> (user_data);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  // The test forest consists of exactly one local tree.
  EXPECT_EQ (ltreeid, 0);
  EXPECT_NE (element, nullptr);
  EXPECT_NE (leaf_elements, nullptr);

  // Count every visited element, including intermediate elements.
  ++(data->callback_count);

  // Verify that the current element face lies on the requested root face.
  EXPECT_TRUE (scheme->element_is_root_boundary (data->eclass, element, face));
  EXPECT_EQ (scheme->element_get_tree_face (data->eclass, element, face), data->root_face);

  if (is_leaf) {
    // Verify that the reported leaf index lies within the local tree.
    EXPECT_GE (tree_leaf_index, 0);
    EXPECT_LT (tree_leaf_index, t8_forest_get_tree_num_leaf_elements (forest, ltreeid));

    ++(data->leaf_count);

    // Inserting the index must succeed because each leaf is visited only once.
    EXPECT_TRUE (data->leaf_indices.insert (tree_leaf_index).second);
  }

  return !data->stop_recursion;
}

/**
 * Parameterized GoogleTest fixture for t8_forest_iterate_faces.
 *
 * Each test configuration specifies an element class, a uniform
 * refinement level, a root face and a callback termination behavior.
 *
 * The fixture constructs a single-tree coarse mesh and a uniformly
 * refined forest using sc_MPI_COMM_SELF. Consequently, the traversal
 * is tested locally and does not require communication between MPI ranks.
 *
 * A separate root element is allocated using the forest's refinement
 * scheme. The root element serves as the starting point for the
 * face iteration.
 *
 * The test fixture also provides a reference implementation that
 * identifies all forest leaves touching a specified root face.
 *
 * The four parameters are:
 * - The element class of the coarse mesh tree.
 * - The uniform refinement level of the forest.
 * - The root face to be traversed.
 * - Whether the callback should terminate recursion immediately.
 */
struct t8_gtest_forest_iterate_faces: public ::testing::TestWithParam<std::tuple<t8_eclass_t, int, int, bool>>
{
 protected:
  /** The uniformly refined forest used by the test. */
  t8_forest_t forest = nullptr;

  /** The refinement scheme associated with the forest. */
  const t8_scheme *scheme = nullptr;

  /** The allocated root element from which face iteration begins. */
  t8_element_t *root = nullptr;

  /** The element class of the forest's single coarse mesh tree. */
  t8_eclass_t eclass = T8_ECLASS_INVALID;

  /**
   * Construct the forest and root element for the current test case.
   *
   * The coarse mesh contains exactly one tree of the element class
   * specified by the test parameters. A uniformly refined forest
   * is then constructed at the requested refinement level.
   *
   *
   * Finally, the root element is allocated and initialized through
   * the forest's refinement scheme. The requested face index is
   * checked against the number of faces of the root element.
   */
  void
  SetUp () override
  {
    const t8_eclass_t test_eclass = std::get<0> (GetParam ());
    const int level = std::get<1> (GetParam ());
    const int face = std::get<2> (GetParam ());

    eclass = test_eclass;

    // Skip configurations with a face index outside the element class range.
    if (face >= t8_eclass_num_faces[eclass]) {
      GTEST_SKIP () << "Face " << face << " is invalid for element class " << t8_eclass_to_string[eclass];
    }

    // Construct a coarse mesh containing a single tree of the requested class.
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_class (cmesh, eclass, sc_MPI_COMM_SELF);

    // Uniformly refine the tree without creating ghost elements.
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_SELF);

    ASSERT_NE (forest, nullptr);
    ASSERT_EQ (t8_forest_get_num_local_trees (forest), 1);

    scheme = t8_forest_get_scheme (forest);

    // Allocate an element and initialize it as the root of the tree.
    scheme->element_new (eclass, 1, &root);
    scheme->set_to_root (eclass, root);

    // Verify that the requested root face is valid for this element class.
    ASSERT_GE (face, 0);
    ASSERT_LT (face, scheme->element_get_num_faces (eclass, root));
  }

  /**
   * Release all resources allocated by SetUp.
   *
   * The root element is destroyed using the same refinement scheme
   * through which it was allocated. The forest reference is then
   * released.
   *
   * Null pointer checks ensure that cleanup is also possible if
   * SetUp terminates before all resources have been initialized.
   */
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
   * Determine the tree-local indices of all leaves touching a root face.
   *
   * This function independently constructs a reference set for the
   * face iteration test. It traverses all leaf elements of the tree
   * and examines their faces using the forest's refinement scheme.
   *
   * A leaf is included if at least one of its faces lies on the
   * root boundary and corresponds to the specified root face.
   *
   * Each qualifying leaf index is inserted into a set, ensuring
   * that a leaf is represented only once even if multiple faces
   * satisfy the boundary condition.
   *
   * The resulting set is used to verify that the face iterator
   * visits all and only the expected leaf elements.
   *
   * \param [in] root_face  The index of the root face for which the
   *                        touching leaf elements are collected.
   *
   * \return               A set containing the tree-local indices
   *                       of all leaves touching the specified face.
   */
  std::set<t8_locidx_t>
  get_expected_leaf_indices (const int root_face) const
  {
    std::set<t8_locidx_t> expected;

    const auto *leaves = t8_forest_tree_get_leaf_elements (forest, 0);
    const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest, 0);

    // Examine every leaf of the local tree independently of the face iterator.
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ++ileaf) {
      const t8_element_t *leaf = t8_element_array_index_locidx (leaves, ileaf);
      const int num_faces = scheme->element_get_num_faces (eclass, leaf);

      // Check whether any face of the current leaf belongs to the root face.
      for (int iface = 0; iface < num_faces; ++iface) {
        if (scheme->element_is_root_boundary (eclass, leaf, iface)
            && scheme->element_get_tree_face (eclass, leaf, iface) == root_face) {
          expected.insert (ileaf);

          // The leaf has already been identified as touching the root face.
          break;
        }
      }
    }

    return expected;
  }
};

/**
 * Verify face iteration for a parameterized forest configuration.
 *
 * The test invokes t8_forest_iterate_faces on the root element
 * and records all callback invocations.
 *
 * Two callback behaviors are evaluated:
 *
 * 1. Complete traversal:
 *    The callback always permits recursion. Every leaf touching
 *    the selected root face must be visited exactly once, and
 *    no leaf outside that face may be included.
 *
 * 2. Early termination:
 *    The callback returns zero upon its first invocation.
 *    The iterator must therefore visit only the root element
 *    and must not descend into any children.
 *
 * At refinement level zero, the root element is itself a leaf.
 * Consequently, a terminating callback still records one leaf.
 *
 * For refined forests, a complete traversal must visit at least
 * one intermediate element in addition to the leaf elements.
 */
TEST_P (t8_gtest_forest_iterate_faces, iterate_faces)
{
  const auto [test_eclass, level, face, stop_recursion] = GetParam ();

  const auto *leaves = t8_forest_tree_get_leaf_elements (forest, 0);

  // Initialize the callback configuration and reset all traversal counters.
  t8_face_test_data data { test_eclass, face, stop_recursion, 0, 0, {} };

  // Determine the expected leaf indices without using the face iterator.
  const auto expected_indices = get_expected_leaf_indices (face);

  // Traverse the requested root face and collect the callback results.
  t8_forest_iterate_faces (forest, 0, root, face, leaves, 0, t8_test_face_callback, &data);

  if (stop_recursion) {
    // A terminating callback must prevent recursion beyond the root element.
    EXPECT_EQ (data.callback_count, 1);

    // At level zero, the root element is already a leaf.
    EXPECT_EQ (data.leaf_count, level == 0 ? 1 : 0);
  }
  else {
    // The iterator must report exactly the expected face leaves.
    EXPECT_EQ (data.leaf_indices, expected_indices);
    EXPECT_EQ (data.leaf_count, expected_indices.size ());

    if (level > 0) {
      // A refined tree requires traversal through intermediate elements.
      EXPECT_GT (data.callback_count, data.leaf_count);
    }
    else {
      // Without refinement, the root is the only visited element.
      EXPECT_EQ (data.callback_count, 1);
    }
  }
}

/**
 * Instantiate the face iteration tests for different element classes,
 * refinement levels, root faces and callback termination behaviors.
 * 
 * Test configurations for:
 * - Element classes: LINE, QUAD and HEX.
 * - Refinement levels: 0 through 3.
 * - Root faces: 0 through 5.
 * - Callback behavior: Complete traversal or immediate termination.
 *
 * Since the number of valid faces depends on the element class,
 * configurations with invalid face indices are skipped in SetUp.
 */
INSTANTIATE_TEST_SUITE_P (AllFaceConfigurations, t8_gtest_forest_iterate_faces,
                          ::testing::Combine (::testing::Values (T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_HEX),
                                              ::testing::Range (0, 4), ::testing::Range (0, 6), ::testing::Bool ()));
