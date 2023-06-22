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
 * test yet. 
 * 
 */
TEST (t8_cmesh_vtk_reader, dummy_test)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh =
    t8_cmesh_vtk_reader ("non-existing-file.vtu", 0, 0, sc_MPI_COMM_WORLD,
                         VTK_FILE_ERROR);
  EXPECT_TRUE (cmesh == NULL);

  cmesh =
    t8_cmesh_vtk_reader ("non-existing-file.vtu", 0, 0, sc_MPI_COMM_WORLD,
                         VTK_UNSTRUCTURED_FILE);
  EXPECT_TRUE (cmesh == NULL);

  cmesh =
    t8_cmesh_vtk_reader ("non-existing-file.vtp", 0, 0, sc_MPI_COMM_WORLD,
                         VTK_POLYDATA_FILE);
  EXPECT_TRUE (cmesh == NULL);

  cmesh =
    t8_cmesh_vtk_reader ("test/testfiles/test_vtk_tri.vtu", 0, 0,
                         sc_MPI_COMM_WORLD, VTK_UNSTRUCTURED_FILE);
  EXPECT_FALSE (cmesh == NULL);
  t8_cmesh_destroy (&cmesh);
  cmesh =
    t8_cmesh_vtk_reader ("test/testfiles/test_vtk_cube.vtp", 0, 0,
                         sc_MPI_COMM_WORLD, VTK_POLYDATA_FILE);
  EXPECT_FALSE (cmesh == NULL);
  t8_cmesh_destroy (&cmesh);

#else
#endif
}
