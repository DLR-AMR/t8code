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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <test/t8_gtest_macros.hxx>

#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"

/* We create and commit a cmesh, then derive a new cmesh
 * from it without any changes.
 * We test whether the new and original cmesh are equal.
 */

/* Note: This test currently fails on many cmeshes and is thus deavtivated.
 * See: https://github.com/DLR-AMR/t8code/issues/920
 */

/* Remove `DISABLED_` from the name of the Test(suite) or use `--gtest_also_run_disabled_tests` when you start working on the issue. */
class DISABLED_t8_cmesh_copy: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    /* Skip test since cmesh copy is not yet working. See https://github.com/DLR-AMR/t8code/issues/920 */

    cmesh_original = GetParam ()->cmesh_create ();

    /* Initialized test cmesh that we derive in the test */
    t8_cmesh_init (&cmesh);
  }

  void
  TearDown () override
  {
    /* Skip test since cmesh copy is not yet working. See https://github.com/DLR-AMR/t8code/issues/920 */

    /* Unref both cmeshes */
    t8_cmesh_unref (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_cmesh_t cmesh_original;
};

static void
test_cmesh_committed (t8_cmesh_t cmesh)
{
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed.";
  ASSERT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) << "Cmesh face consistency failed.";
}

TEST_P (DISABLED_t8_cmesh_copy, test_cmesh_copy)
{
  t8_cmesh_set_derive (cmesh, cmesh_original);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  test_cmesh_committed (cmesh);

  EXPECT_TRUE (t8_cmesh_is_equal (cmesh, cmesh_original));
}

/* Test all cmeshes over all different inputs*/

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_copy, DISABLED_t8_cmesh_copy, AllCmeshsParam, pretty_print_base_example);
