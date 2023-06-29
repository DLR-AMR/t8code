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

/**
 * This is currently a place-holder for a propper cmesh_vtk_reader-test.
 * The function is not implemented yet and therefore we do not provide a proper
 * test yet. A proper test would compare the read files with a reference-cmesh. 
 * 
 */


/* *INDENT-OFF* */
class vtk_reader : public testing::TestWithParam<vtk_file_type_t>{
  protected:
    void SetUp() override{
      file_type = GetParam();
      file = (int)file_type + 1;
    }
    int file;
    vtk_file_type_t file_type;
    const char* failing_files[3] = {
      "no_file",
      "non-existing-file.vtu",
      "non-existing-file.vtp"
    };
    const char* test_files[3] = {
      "no_file",
      "test/testfiles/test_vtk_tri.vtu",
      "test/testfiles/test_vtk_cube.vtp"
    };
    const int num_points[3] = {0, 121, 24};
};
/* *INDENT-ON* */

/* All readers should fail properly with a non-existing file. */
TEST_P (vtk_reader, vtk_to_cmesh_fail)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh =
    t8_cmesh_vtk_reader (failing_files[file], 0, 0, sc_MPI_COMM_WORLD,
                         file_type);
  EXPECT_TRUE (cmesh == NULL);
#else
#endif
}

/* All readers should construct a cmesh from a file. */
TEST_P (vtk_reader, vtk_to_cmesh_success)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh = t8_cmesh_vtk_reader (test_files[file], 0, 0,
                                                   sc_MPI_COMM_WORLD,
                                                   file_type);
  if (file_type != VTK_FILE_ERROR) {
    EXPECT_FALSE (cmesh == NULL);
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
    vtkSmartPointer < vtkPointSet > points =
      t8_vtk_reader_pointSet (test_files[file], 0, 0,
                              sc_MPI_COMM_WORLD, file_type);
    int                 test_points = points->GetNumberOfPoints ();
    EXPECT_EQ (num_points[file], test_points);
  }
#else
#endif
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_vtk_reader, vtk_reader,
                          testing::Range (VTK_FILE_ERROR, VTK_NUM_TYPES));
