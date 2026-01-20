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
#include <t8_cmesh/t8_cmesh.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_partition.h>
#include <test/t8_gtest_macros.hxx>

#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

/* We create and commit a cmesh, then derive a new cmesh
 * from it without any changes.
 * We test whether the new and original cmesh are equal.
 */

class t8_cmesh_copy: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh_original = GetParam ()->cmesh_create ();
  }

  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh_original);
  }
  t8_cmesh_t cmesh_original;
};

static void
test_cmesh_committed (t8_cmesh_t cmesh)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed.";
  ASSERT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) << "Cmesh face consistency failed.";
}

TEST_P (t8_cmesh_copy, test_cmesh_copy)
{
  t8_cmesh_t cmesh_copy;
  t8_cmesh_init (&cmesh_copy);
  t8_cmesh_ref (cmesh_original);
  t8_cmesh_set_derive (cmesh_copy, cmesh_original);
  t8_cmesh_commit (cmesh_copy, sc_MPI_COMM_WORLD);

  test_cmesh_committed (cmesh_copy);
  EXPECT_TRUE (t8_cmesh_is_equal (cmesh_copy, cmesh_original));
  t8_cmesh_unref (&cmesh_copy);
}

/* Test all cmeshes over all different inputs*/
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_copy, t8_cmesh_copy, AllCmeshsParam, pretty_print_base_example);
