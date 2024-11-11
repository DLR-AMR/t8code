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
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

class forest_search: public testing::TestWithParam<std::tuple<t8_eclass, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());

    default_scheme = t8_scheme_new_default ();
    /* Construct a cube coarse mesh */
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    /* Build a uniform forest */
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
  int level;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
  t8_scheme *default_scheme;
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
  if (is_leaf) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, tree_class);
    const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    /* Set the corresponding entry to 1 */
    (*user_data)[tree_offset + tree_leaf_index] = true;
    /* Test whether tree_leaf_index is actually the index of the element */
    t8_locidx_t test_ltreeid;
    const t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);

    EXPECT_ELEM_EQ (scheme, tree_class, element, test_element);
    EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
  }
  return true;
}

void
t8_test_search_query_all_fn (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                             const bool is_leaf, const t8_element_array_t *leaf_elements,
                             const t8_locidx_t tree_leaf_index, std::vector<int> &queries,
                             std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches,
                             std::vector<bool> *user_data)
{
  EXPECT_FALSE (queries.empty ()) << "query callback must be called with queries argument. ";
  EXPECT_EQ (active_query_indices.size (), (long unsigned int) 1)
    << "Wrong number of active queries passed to query callback.";
  for (int iquery : active_query_indices) {
    /* The query is an int with value 42 (see below) */
    EXPECT_EQ (queries[iquery], 42) << "Wrong query argument passed to query callback.";
    if (is_leaf) {
      /* Test whether tree_leaf_index is actually the index of the element */
      t8_locidx_t test_ltreeid;
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
      const t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, tree_class);

      const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
      const t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);
      EXPECT_ELEM_EQ (ts, element, test_element);
      EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
    }
    query_matches[iquery] = true;
  }
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
  std::for_each (matched_leaves.begin (), matched_leaves.end (),
                 [] (bool b) { ASSERT_TRUE (b) << "Search did not match all leaves. First mismatch at leaf " << b; });

  t8_forest_unref (&forest);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_search, forest_search, testing::Combine (AllEclasses, testing::Range (0, 6)));
