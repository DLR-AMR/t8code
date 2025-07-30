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
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include "t8_cmesh/t8_cmesh_trees.h"

/* In this file we test the msh file (gmsh) reader of the cmesh.
 * Currently, we support version 2 and 4 ascii.
 * We read a mesh file and check whether the constructed cmesh is correct.
 * We also try to read version 2 binary and version 4 binary
 * formats. All are not supported and we expect the reader to catch this.
 */

static void
t8_supported_msh_file (t8_cmesh_t cmesh)
{
  double *vertices;
  t8_locidx_t ltree_id;
  t8_locidx_t lnum_trees;
  t8_eclass_t tree_class;

  /* Description of the properties of the example msh-files. */
  const int number_elements = 4;
  const t8_eclass_t elem_type = T8_ECLASS_TRIANGLE;

  int vertex[6][2] = { { 0, 0 }, { 2, 0 }, { 4, 0 }, { 1, 2 }, { 3, 2 }, { 2, 4 } };

  int elements[4][3] = { { 0, 1, 3 }, { 1, 4, 3 }, { 1, 2, 4 }, { 3, 4, 5 } };

  int face_neigh_elem[4][3] = { { 1, -1, -1 }, { 3, 0, 2 }, { -1, 1, -1 }, { -1, -1, 1 } };

  ASSERT_FALSE (cmesh == NULL) << "Reading cmesh failed.";

  /* Checks if the cmesh was committed. */
  ASSERT_TRUE (t8_cmesh_is_committed (cmesh)) << "Cmesh commit failed";
  /* Checks for face consistency. */
  ASSERT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees)) << "Cmesh face consistency failed.";

  /* Checks if the number of elements was read correctly. */
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), number_elements) << "Number of elements in msh-file was read incorrectly.";

  /* Number of local trees. */
  lnum_trees = t8_cmesh_get_num_local_trees (cmesh);
  /* Iterate through the local elements and check if they were read properly. */
  /* All trees should be local to the master rank. */
  for (t8_locidx_t ltree_it = 0; ltree_it < lnum_trees; ltree_it++) {
    tree_class = t8_cmesh_get_tree_class (cmesh, ltree_it);
    ASSERT_FALSE (t8_eclass_compare (tree_class, elem_type)) << "Element type in msh-file was read incorrectly.";

    /* Get pointer to the vertices of the tree. */
    vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_it);
    /* Checking the msh-files elements and nodes. */
    for (int i = 0; i < 3; i++) {
      /* Checks if x and y coordinate of the nodes are not read correctly. */
      ASSERT_EQ (vertex[elements[ltree_it][i]][0], (int) vertices[3 * i]) << "x coordinate was read incorrectly";
      ASSERT_EQ (vertex[elements[ltree_it][i]][1], (int) vertices[(3 * i) + 1]) << "y coordinate was read incorrectly";

      /* Checks whether the face neighbor elements are not read correctly. */
      ltree_id = t8_cmesh_get_face_neighbor (cmesh, ltree_it, i, NULL, NULL);
      ASSERT_EQ (ltree_id, face_neigh_elem[ltree_it][i])
        << "The face neighbor element in the example test file was not read correctly.";
      const t8_eclass_t neighbor_eclass = t8_cmesh_get_tree_face_neighbor_eclass (cmesh, ltree_it, i);
      // If a face neighbor exists, the return value must match the element type, otherwise it must
      // be T8_ECLASS_INVALID
      const t8_eclass_t reference_value = face_neigh_elem[ltree_it][i] == -1 ? T8_ECLASS_INVALID : elem_type;
      EXPECT_EQ (neighbor_eclass, reference_value) << "mismatch in face neighbor eclass.";
    }
  }
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers2_ascii)
{

  const char fileprefix[BUFSIZ - 4] = "test/testfiles/test_msh_file_vers2_ascii";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_debugf ("Checking msh file version 2 ascii...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;

  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from ascii version 2, but should be able to.";

  t8_supported_msh_file (cmesh);

  /* The cmesh was read successfully and we need to destroy it. */
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_ascii)
{

  const char fileprefix[BUFSIZ - 4] = "test/testfiles/test_msh_file_vers4_ascii";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_debugf ("Checking msh file version 4 ascii...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from ascii version 4, but should be able to.";

  /* The cmesh was read successfully and we need to destroy it. */
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers2_bin)
{

  const char fileprefix[BUFSIZ - 4] = "test/testfiles/test_msh_file_vers2_bin";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_debugf ("Checking msh file version 2 binary...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh == NULL) << "Expected fail of reading binary msh file v.2, but did not fail.";

  t8_debugf ("Error handling successful.\n");
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_bin)
{

  const char fileprefix[BUFSIZ - 4] = "test/testfiles/test_msh_file_vers4_bin";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_debugf ("Checking msh file version 4 binary...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh == NULL) << "Expected fail of reading binary msh file v.4, but did not fail.";

  t8_debugf ("Error handling successful.\n");
}
