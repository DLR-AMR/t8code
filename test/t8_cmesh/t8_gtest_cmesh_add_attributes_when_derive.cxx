/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include "t8_cmesh_generator/t8_cmesh_example_sets.hxx"

/* Test if adding attributes to a cmesh when deriving it from another
 * cmesh is working.
 * This test is currently disabled, see https://github.com/DLR-AMR/t8code/issues/923 */

class DISABLED_t8_cmesh_add_attributes: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh = GetParam ()->cmesh_create ();

    t8_cmesh_init (&cmesh_derived);
    t8_cmesh_set_derive (cmesh_derived, cmesh);

    const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
    const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh);
    const t8_locidx_t num_local_and_ghost = num_local_trees + num_ghost_trees;

    /* For each local tree and ghost, we set two attributes. */
    for (t8_locidx_t itree = 0; itree < num_local_and_ghost; ++itree) {
      t8_locidx_t locidx_attribute = itree;
      const t8_gloidx_t gtreeid = t8_cmesh_get_global_id (cmesh, itree);
      const int data_persists = 0;

      t8_cmesh_set_attribute (cmesh_derived, gtreeid, t8_testsuite_get_package_id (), 0, &locidx_attribute,
                              sizeof (t8_locidx_t), data_persists);
      t8_cmesh_set_attribute_string (cmesh_derived, gtreeid, t8_testsuite_get_package_id (), 1, string_attribute);
    }
    t8_cmesh_commit (cmesh_derived, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh_derived);
  }

  t8_cmesh_t cmesh;
  t8_cmesh_t cmesh_derived;
  const char *string_attribute = "This is a test attribute";
};

/** Check attribute values of cmeshes against reference values. */
TEST_P (DISABLED_t8_cmesh_add_attributes, check_attributes)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh_derived);
  const t8_locidx_t num_ghost_trees = t8_cmesh_get_num_ghosts (cmesh_derived);
  const t8_locidx_t num_local_and_ghost = num_local_trees + num_ghost_trees;

  for (t8_locidx_t itree = 0; itree < num_local_and_ghost; ++itree) {
    const t8_locidx_t check_locidx_attribute
      = *(t8_locidx_t *) t8_cmesh_get_attribute (cmesh_derived, t8_testsuite_get_package_id (), 0, itree);

    const char *check_string_attribute
      = (const char *) t8_cmesh_get_attribute (cmesh_derived, t8_testsuite_get_package_id (), 1, itree);
    EXPECT_EQ (itree, check_locidx_attribute);
    EXPECT_STREQ (check_string_attribute, string_attribute);
  }
}

/* Test for different number of trees. */
INSTANTIATE_TEST_SUITE_P (t8_gtest_add_attributes_when_derive, DISABLED_t8_cmesh_add_attributes, AllCmeshsParam,
                          pretty_print_base_example);
