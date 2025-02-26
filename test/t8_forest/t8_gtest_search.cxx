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
#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_search/t8_forest_search.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_schemes.hxx>

class forest_search: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass>, int>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    level = std::get<1> (GetParam ());

    /* Construct a cube coarse mesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    /* Build a uniform forest */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
  int level;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  const t8_scheme *scheme;
};

/* A search function that matches all elements.
 * This function assumes that the forest user pointer is an sc_array
 * with one int for each local leaf.
 * If this function is called for a leaf, it sets the corresponding entry to 1.
 */
bool
t8_test_search_all_fn (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                       const bool is_leaf, const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index,
                       std::vector<bool> *user_data)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (user_data != nullptr);
  if (is_leaf) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    const t8_scheme *ts = t8_forest_get_scheme (forest);
    const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    /* Set the corresponding entry to 1 */
    (*user_data)[tree_offset + tree_leaf_index] = true;
    /* Test whether tree_leaf_index is actually the index of the element */
    t8_locidx_t test_ltreeid;
    const t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);

    EXPECT_ELEM_EQ (ts, tree_class, element, test_element);
    EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
  }
  return true;
}

inline bool
t8_test_search_query_all_fn (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                             const bool is_leaf, const t8_element_array_t *leaf_elements,
                             const t8_locidx_t tree_leaf_index, const int &querie, std::vector<bool> *user_data)
{
  /* The query is an int with value 42 (see below) */
  EXPECT_EQ (querie, 42) << "Wrong query argument passed to query callback.";
  if (is_leaf) {
    /* Test whether tree_leaf_index is actually the index of the element */
    t8_locidx_t test_ltreeid;
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    const t8_scheme *ts = t8_forest_get_scheme (forest);

    const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    const t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);
    EXPECT_ELEM_EQ (ts, tree_class, element, test_element);
    EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
  }
  return true;
}

TEST_P (forest_search, t8_test_search_all_fn)
{
  t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  /* set up an array in which we flag whether an element was matched in the
   * search */
  std::vector<bool> matched_leaves (num_elements, false);

  /* Call search. This search matches all elements. After this call we expect
   * all entries in the matched_leaves array to be set to 1. */

  t8_search<std::vector<bool>> search (t8_test_search_all_fn);

  search.update_user_data (&matched_leaves);
  search.update_forest (forest);
  search.do_search ();

  /* Check whether matched_leaves entries are all 1 */
  for (size_t i = 0; i < matched_leaves.size (); ++i) {
    ASSERT_TRUE (matched_leaves[i]) << "Search did not match all leaves. Mismatch at leaf " << i;
  }

  t8_forest_unref (&forest);
}

TEST_P (forest_search, test_search_one_query_matches_all)
{
  /* set up a single query containing our query */
  std::vector<int> queries = { 42 };

  t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  /* set up an array in which we flag whether an element was matched in the
   * search */
  std::vector<bool> matched_leaves (num_elements, false);

  /* Call search. This search matches all elements. After this call we expect
   * all entries in the matched_leaves array to be set to 1. */

  t8_search_with_queries<int, std::vector<bool>> search (t8_test_search_all_fn, t8_test_search_query_all_fn, queries);

  search.update_queries (queries);
  search.update_user_data (&matched_leaves);
  search.update_forest (forest);
  search.do_search ();

  /* Check whether matched_leaves entries are all 1 */
  for (size_t i = 0; i < matched_leaves.size (); ++i) {
    ASSERT_TRUE (matched_leaves[i]) << "Search did not match all leaves. Mismatch at leaf " << i;
  }

  t8_forest_unref (&forest);
}
#if T8CODE_TEST_LEVEL >= 2
const int maxlvl = 5;
#else
const int maxlvl = 6;
#endif

INSTANTIATE_TEST_SUITE_P (t8_gtest_search, forest_search, testing::Combine (AllSchemes, testing::Range (0, maxlvl)));
