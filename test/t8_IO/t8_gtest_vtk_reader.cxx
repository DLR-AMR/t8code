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
#include <t8_cmesh_vtk_reader.hxx>
#include <test/t8_gtest_macros.hxx>

#define T8_VTK_TEST_NUM_PROCS 2

/**
 * This is currently a place-holder for a proper cmesh_vtk_reader-test.
 * The function is not implemented yet and therefore we do not provide a proper
 * test yet. A proper test would compare the read files with a reference-cmesh. 
 * 
 */
const vtk_file_type_t gtest_vtk_filetypes[VTK_NUM_TYPES]
  = { VTK_FILE_ERROR, VTK_UNSTRUCTURED_FILE, VTK_POLYDATA_FILE, VTK_PARALLEL_UNSTRUCTURED_FILE,
      VTK_PARALLEL_POLYDATA_FILE };

class vtk_reader: public testing::TestWithParam<std::tuple<int, int, int>> {
 protected:
  void
  SetUp () override
  {
    int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);

    /* `mpisize` must match with the number of files that are read in. */
    if (mpisize != T8_VTK_TEST_NUM_PROCS) {
      GTEST_SKIP ();
    }
    file = std::get<0> (GetParam ());
    file_type = gtest_vtk_filetypes[file];
    partition = std::get<1> (GetParam ());
    main_proc = std::get<2> (GetParam ());
    distributed = (file_type & VTK_PARALLEL_FILE) && partition;
  }
  int file;
  vtk_file_type_t file_type;
  int partition;
  int distributed;
  int main_proc;
  int mpisize;
  int mpirank;
  const char* failing_files[5] = { "no_file", "non-existing-file.vtu", "non-existing-file.vtp",
                                   "non-existing-file.pvtu", "non-existing-file.pvtp" };
  const char* test_files[5] = { "no_file", "test/testfiles/test_vtk_tri.vtu", "test/testfiles/test_vtk_cube.vtp",
                                "test/testfiles/test_parallel_file.pvtu", "test/testfiles/test_polydata.pvtp" };
  const int num_points[5] = { 0, 121, 24, 6144, 900 };
  const int num_trees[5] = { 0, 200, 12, 1024, 1680 };
};

/* All readers should fail properly with a non-existing file. */
TEST_P (vtk_reader, vtk_to_cmesh_fail)
{
#if T8_WITH_VTK
  t8_cmesh_t cmesh
    = t8_cmesh_vtk_reader (failing_files[file], 0, main_proc, sc_MPI_COMM_WORLD, file_type, t8_testsuite_package_id, 0);
  EXPECT_TRUE (cmesh == NULL);
#else
#endif
}

/* All readers should construct a cmesh from a file. */
TEST_P (vtk_reader, vtk_to_cmesh_success)
{
#if T8_WITH_VTK
  int mpirank;
  int mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
  t8_cmesh_t cmesh = t8_cmesh_vtk_reader (test_files[file], partition, main_proc, sc_MPI_COMM_WORLD, file_type,
                                          t8_testsuite_package_id, 0);
  if (file_type != VTK_FILE_ERROR) {
    EXPECT_FALSE (cmesh == NULL);
    const int test_num_trees = t8_cmesh_get_num_local_trees (cmesh);
    if (distributed) {
      /* In this testcase the cmesh should be distributed equally. */
      EXPECT_EQ (num_trees[file] / mpisize, test_num_trees);
    }
    else {
      if (!partition || main_proc == mpirank) {
        /* The proc has the complete cmesh */
        EXPECT_EQ (num_trees[file], test_num_trees);
      }
      else {
        /* Every other proc should be empty. */
        EXPECT_EQ (0, test_num_trees);
      }
    }
    t8_cmesh_destroy (&cmesh);
  }
  else {
    EXPECT_TRUE (cmesh == NULL);
  }
#else
#endif
}

/* Read a file as a pointSet and compare the number of points with the known number of points. */
TEST_P (vtk_reader, vtk_to_pointSet)
{
#if T8_WITH_VTK
  if (file_type != VTK_FILE_ERROR) {
    vtkSmartPointer<vtkPointSet> points
      = t8_vtk_reader_pointSet (test_files[file], partition, main_proc, sc_MPI_COMM_WORLD, file_type);
    int test_points = points->GetNumberOfPoints ();
    if (distributed) {
      /* The points should be distributed equally in this case. */
      EXPECT_EQ (num_points[file] / mpisize, test_points);
    }
    else {
      if (!partition || main_proc == mpirank) {
        /* The proc has all points. */
        EXPECT_EQ (num_points[file], test_points);
      }
      else {
        /* Every other proc should have no points. */
        EXPECT_EQ (0, test_points);
      }
    }
  }
#else
#endif
}

/* Currently does not work for parallel files. Replace with VTK_NUM_TYPES as soon
 * as reading and constructing cmeshes from parallel files is enabled. */
INSTANTIATE_TEST_SUITE_P (t8_gtest_vtk_reader, vtk_reader,
                          testing::Combine (testing::Range (VTK_FILE_ERROR + 1, (int) VTK_NUM_TYPES),
                                            testing::Values (0, 1), testing::Range (0, T8_VTK_TEST_NUM_PROCS)));
