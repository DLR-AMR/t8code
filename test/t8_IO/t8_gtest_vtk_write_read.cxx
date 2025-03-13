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

class vtk_write_read_test: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    cmesh = GetParam ()->cmesh_create (); /* create a cmesh from example base */
  }

  void
  TearDown () override
  {
    std::remove ("test_vtk_0.vtu");
    std::remove ("test_vtk.pvtu"); /* delete files after test finished */
  }
  t8_cmesh_t cmesh;
};

TEST_P (vtk_write_read_test, write_read_vtk)
{
  t8_cmesh_t unpartitioned_cmesh_read, unpartitioned_cmesh_write;

  const int write_success = t8_cmesh_vtk_write_file_via_API (cmesh, "test_vtk", sc_MPI_COMM_WORLD);
  const bool partitioned = t8_cmesh_is_partitioned (cmesh);
  t8_cmesh_t cmesh2 = t8_vtk_reader_cmesh ("test_vtk.pvtu", partitioned, 0, sc_MPI_COMM_WORLD, VTK_SERIAL_FILE);
  ASSERT_EQ (cmesh2 != NULL, write_success);

  if (cmesh2 != NULL) {
    if (partitioned) { /* repartition if cmesh is partitioned */
      t8_cmesh_init (&unpartitioned_cmesh_write);
      t8_cmesh_set_derive (unpartitioned_cmesh_write, cmesh);
      t8_cmesh_set_partition_offsets (unpartitioned_cmesh_write,
                                      t8_cmesh_offset_percent (cmesh, sc_MPI_COMM_WORLD, 100));
      t8_cmesh_commit (unpartitioned_cmesh_write, sc_MPI_COMM_WORLD);
      t8_cmesh_init (&unpartitioned_cmesh_read);
      t8_cmesh_set_derive (unpartitioned_cmesh_read, cmesh2);
      t8_cmesh_set_partition_offsets (unpartitioned_cmesh_read,
                                      t8_cmesh_offset_percent (cmesh2, sc_MPI_COMM_WORLD, 100));
      t8_cmesh_commit (unpartitioned_cmesh_read, sc_MPI_COMM_WORLD);

      /* check equality of repartioned cmeshes */
      t8_gloidx_t num_trees_cmesh_1 = t8_cmesh_get_num_trees (unpartitioned_cmesh_write);
      t8_gloidx_t num_trees_cmesh_2 = t8_cmesh_get_num_trees (unpartitioned_cmesh_read);
      ASSERT_EQ (num_trees_cmesh_1, num_trees_cmesh_2);
      EXPECT_TRUE (t8_cmesh_is_equal_ext (unpartitioned_cmesh_write, unpartitioned_cmesh_read, 0));

      /* destroy both cmeshes */
      t8_cmesh_unref (&unpartitioned_cmesh_read);
      t8_cmesh_unref (&unpartitioned_cmesh_write);
    }
    else { /* continue normally without repartitioning */
      t8_gloidx_t num_trees_cmesh_1 = t8_cmesh_get_num_trees (cmesh);
      t8_gloidx_t num_trees_cmesh_2 = t8_cmesh_get_num_trees (cmesh2);
      ASSERT_EQ (num_trees_cmesh_1, num_trees_cmesh_2);
      EXPECT_TRUE (t8_cmesh_is_equal_ext (cmesh, cmesh2, 0));
      t8_cmesh_unref (&cmesh2);
      t8_cmesh_unref (&cmesh);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (Test_vtk_write_read, vtk_write_read_test, AllCmeshsParam, pretty_print_base_example);
