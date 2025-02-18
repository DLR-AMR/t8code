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
#include <t8_vtk/t8_vtk_writer.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <t8_cmesh/t8_cmesh_partition.h>

#include <t8_vtk/t8_vtk_writer.h>
#include <t8_vtk/t8_vtk_reader.hxx>

#define EXPECT_CMESH_EQ(cmesh_a, cmesh_b) EXPECT_TRUE (t8_cmesh_is_equal_ext (cmesh_a, cmesh_b, 0))

class vtk_write_read_test: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh_write = GetParam ()->cmesh_create (); /* create a cmesh from example base */
  }

  void
  TearDown () override
  {
    std::remove ("test_vtk_0.vtu");
    std::remove ("test_vtk.pvtu"); /* delete files after test finished */
    t8_cmesh_destroy (&cmesh_write);
  }
  t8_cmesh_t cmesh_write;
};

TEST_P (vtk_write_read_test, write_read_vtk)
{
  const int write_success = t8_cmesh_vtk_write_file_via_API (cmesh_write, "test_vtk", sc_MPI_COMM_WORLD);
  const int partitioned = t8_cmesh_is_partitioned (cmesh_write);
  t8_cmesh_t cmesh_read = t8_vtk_reader_cmesh ("test_vtk.pvtu", partitioned, 0, sc_MPI_COMM_WORLD, VTK_SERIAL_FILE);
  ASSERT_EQ (cmesh_read != NULL, write_success);

  if (cmesh_read != NULL) {
    /* check equality of cmeshes */
    const t8_gloidx_t num_trees_cmesh_1 = t8_cmesh_get_num_trees (cmesh_write);
    const t8_gloidx_t num_trees_cmesh_2 = t8_cmesh_get_num_trees (cmesh_read);
    ASSERT_EQ (num_trees_cmesh_1, num_trees_cmesh_2);
    EXPECT_CMESH_EQ (cmesh_write, cmesh_read);
    t8_cmesh_destroy (&cmesh_read);
  }
}

INSTANTIATE_TEST_SUITE_P (Test_vtk_write_read, vtk_write_read_test, AllCmeshsParam, pretty_print_base_example);
