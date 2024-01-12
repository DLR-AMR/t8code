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
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>

/* Test if multiple attributes are partitioned correctly. */

#define T8_ATTRIBUTE_TEST_MAX_NUM_ENTRIES 1000

class cmesh_attribute_gloidx_array: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    num_entries = GetParam ();

    /* Build a cmesh with one QUAD tree and one TRIANGLE tree. */
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_dimension (cmesh, 2);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
    t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);

    /* Allocate space for entries. */
    entries = T8_ALLOC (t8_gloidx_t, num_entries);

    /* Fill with 0, 1, 2, 3, 4 ... */
    for (t8_locidx_t ientry = 0; ientry < num_entries; ++ientry) {
      entries[ientry] = ientry;
    }

    t8_debugf ("Calling set_attribute with count %i\n", num_entries);
    /* Set the array as attribute twice. Once with data_persist and once without. */
    /* Attribute at tree 0, data_persist = 0 */
    t8_gloidx_t tree_with_attribute = 0;
    int data_persists = 0;
    t8_cmesh_set_attribute_gloidx_array (cmesh, tree_with_attribute, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY,
                                         entries, num_entries, data_persists);

    /* Attribute at tree 1, data_persist = 1 */
    tree_with_attribute = 1;
    data_persists = 1;
    t8_cmesh_set_attribute_gloidx_array (cmesh, tree_with_attribute, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY,
                                         entries, num_entries, data_persists);

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
    /* It is save to free the entries after commit, since the value got copied. */
    T8_FREE (entries);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_locidx_t num_entries;
  t8_gloidx_t *entries;
  t8_gloidx_t *get_entries;
};

/** Check attribute values of cmeshes against reference values. */
TEST_P (cmesh_attribute_gloidx_array, check_values_data_not_persist)
{
  const t8_locidx_t tree_id = 0;
  get_entries = t8_cmesh_get_attribute_gloidx_array (cmesh, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY, tree_id,
                                                     num_entries);

  /* If we did not store any values, we except to get the NULL pointer back. */
  if (entries == NULL) {
    EXPECT_TRUE (get_entries == NULL);
    /* Number of entries must be < 0 in that case and we cannot continue the test otherwise. */
    ASSERT_EQ (num_entries, 0);
  }
  else {
    /* Otherwise it must not be NULL and we abort the test if it is. */
    ASSERT_TRUE (get_entries != NULL);
  }

  /* Check for equality of the values. */
  for (t8_locidx_t ientry = 0; ientry < num_entries; ++ientry) {
    EXPECT_EQ (get_entries[ientry], ientry);
  }
}

/** Check attribute values of cmeshes against reference values. */
TEST_P (cmesh_attribute_gloidx_array, check_values_data_persist)
{
  const t8_locidx_t tree_id = 1;
  get_entries = t8_cmesh_get_attribute_gloidx_array (cmesh, t8_get_package_id (), T8_CMESH_NEXT_POSSIBLE_KEY, tree_id,
                                                     num_entries);

  /* If we did not store any values, we except to get the NULL pointer back. */
  if (entries == NULL) {
    EXPECT_TRUE (get_entries == NULL);
    /* Number of entries must be < 0 in that case and we cannot continue the test otherwise. */
    ASSERT_EQ (num_entries, 0);
  }
  else {
    /* Otherwise it must not be NULL and we abort the test if it is. */
    ASSERT_TRUE (get_entries != NULL);
  }

  /* Check for equality of the values. */
  for (t8_locidx_t ientry = 0; ientry < num_entries; ++ientry) {
    EXPECT_EQ (get_entries[ientry], ientry);
  }
}

/* Test for different number of entries.
 * 0, 100, 200, ... T8_ATTRIBUTE_TEST_MAX_NUM_ENTRIES */
INSTANTIATE_TEST_SUITE_P (t8_gtest_attribute_gloidx_array, cmesh_attribute_gloidx_array,
                          testing::Range (0, T8_ATTRIBUTE_TEST_MAX_NUM_ENTRIES + 1, 100));
