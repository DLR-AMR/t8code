/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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
#include <unistd.h> /* Needed to check for file access */
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include "t8_test_data_dir.h"
#include <string>

/* In this file we test the msh file (gmsh) reader of the cmesh.
 * Currently, we support version 2 and 4 ascii.
 * We read a mesh file and check whether the constructed cmesh is correct.
 * We also try to read version 2 binary and version 4 binary
 * formats. All are not supported and we expect the reader to catch this.
 */

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_ascii)
{

  std::string fileprefix = std::string (T8_TEST_DATA_DIR) + "/test_msh_file_vers4_ascii";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix.c_str ());

  t8_debugf ("Checking msh file version 4 ascii...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix.c_str (), 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from ascii version 4, but should be able to.";

  /* The cmesh was read successfully and we need to destroy it. */
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_bin)
{

  std::string fileprefix = std::string (T8_TEST_DATA_DIR) + "/test_msh_file_vers4_bin";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix.c_str ());

  t8_debugf ("Checking msh file version 4 binary...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix.c_str (), 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh == NULL) << "Expected fail of reading binary msh file v.4, but did not fail.";

  t8_debugf ("Error handling successful.\n");
}
