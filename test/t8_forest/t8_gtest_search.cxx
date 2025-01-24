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
#include <t8_forest/t8_forest_iterate.h>
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
static int
t8_test_search_all_fn (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                       const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index)
{
  sc_array_t *matched_leaves = (sc_array_t *) t8_forest_get_user_data (forest);
  if (is_leaf) {
    t8_locidx_t test_ltreeid;
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);

    const t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
    /* Set the corresponding entry to 1 */
    *(int *) t8_sc_array_index_locidx (matched_leaves, tree_offset + tree_leaf_index) = 1;
    /* Test whether tree_leaf_index is actually the index of the element */
    const t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);

    EXPECT_ELEM_EQ (scheme, tree_class, element, test_element);
    EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
  }
  return 1;
}

static void
t8_test_search_query_all_fn (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                             const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index,
                             sc_array_t *queries, sc_array_t *query_indices, int *query_matches,
                             const size_t num_active_queries)
{
  EXPECT_TRUE (queries != NULL) << "query callback must be called with queries argument. ";
  EXPECT_EQ (num_active_queries, (long unsigned int) 1) << "Wrong number of active queries passed to query callback.";
  for (size_t iquery = 0; iquery < num_active_queries; iquery++) {
    void *query = sc_array_index_int (queries, iquery);
    /* The query callback is always called with a query */
    EXPECT_TRUE (query != NULL) << "query " << iquery << " is NULL.";
    /* The query is an int with value 42 (see below) */
    EXPECT_EQ (*(int *) query, 42) << "Wrong query argument passed to query callback.";
    if (is_leaf) {
      /* Test whether tree_leaf_index is actually the index of the element */
      t8_locidx_t test_ltreeid;
      const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
      const t8_scheme *scheme = t8_forest_get_scheme (forest);

      t8_locidx_t tree_offset = t8_forest_get_tree_element_offset (forest, ltreeid);
      t8_element_t *test_element = t8_forest_get_element (forest, tree_offset + tree_leaf_index, &test_ltreeid);
      EXPECT_ELEM_EQ (scheme, tree_class, element, test_element);
      EXPECT_EQ (ltreeid, test_ltreeid) << "Tree mismatch in search.";
    }
    query_matches[iquery] = 1;
  }
}

TEST_P (forest_search, test_search_one_query_matches_all)
{
  const int query = 42;
  sc_array_t queries;
  sc_array_t matched_leaves;

  /* set up a single query containing our query */
  sc_array_init_size (&queries, sizeof (int), 1);
  *(int *) sc_array_index (&queries, 0) = query;

  t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);
  /* set up an array in which we flag whether an element was matched in the
   * search */
  sc_array_init_size (&matched_leaves, sizeof (int), num_elements);
  /* write 0 in every entry */
  for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
    *(int *) t8_sc_array_index_locidx (&matched_leaves, ielement) = 0;
  }

  /* Set the array as user data so that we can access it in the search callback */
  t8_forest_set_user_data (forest, &matched_leaves);
  /* Call search. This search matches all elements. After this call we expect
   * all entries in the matched_leaves array to be set to 1. */

  t8_forest_search (forest, t8_test_search_all_fn, t8_test_search_query_all_fn, &queries);

  /* Check whether matched_leaves entries are all 1 */
  for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
    ASSERT_TRUE (*(int *) t8_sc_array_index_locidx (&matched_leaves, ielement))
      << "Search did not match all leaves. First mismatch at leaf " << ielement;
  }

  t8_forest_unref (&forest);
  sc_array_reset (&matched_leaves);
  sc_array_reset (&queries);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_search, forest_search, testing::Combine (AllSchemes, testing::Range (0, 6)));
