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

#include <unistd.h>             /* Needed to check for file access */
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include "t8_cmesh/t8_cmesh_trees.h"

/* In this file we test the msh file (gmsh) reader of the cmesh.
 * Currently, we support version 2 ascii.
 * We read a mesh file and check whether the constructed cmesh is correct.
 * We also try to read version 2 binary, version 4 ascii and version 4 binary
 * formats. All are not supported and we expect the reader to catch this.
 */

/* Check whether the input cmesh matches a given coarse mesh.
 * Returns 1 on success, -1 on failure, 0 if cmesh is NULL.
 * We use this function to check whether the cmesh was read correctly from
 * the file test_msh_file_vers2_ascii.msh
 */
static int
t8_test_supported_msh_file (t8_cmesh_t cmesh)
{
  int                 retval;
  int                 read_node = 1;
  int                 check_neigh_elem = 1;
  int                 check_element_type;
  double             *vertices;
  t8_gloidx_t         num_gtree;
  t8_locidx_t         ltree_id;
  t8_locidx_t         ltree_it;
  t8_locidx_t         lnum_trees;
  t8_eclass_t         tree_class;

  /* Description of the properties of the example msh-files. */
  const int           number_elements = 4;
  const t8_eclass_t   elem_type = T8_ECLASS_TRIANGLE;
  /* *INDENT-OFF* */
  int vertex[6][2] = {
                       {0, 0},
                       {2, 0},
                       {4, 0},
                       {1, 2},
                       {3, 2},
                       {2, 4} };

  int elements[4][3] = {
                        {0, 1, 3},
                        {1, 4, 3},
                        {1, 2, 4},
                        {3, 4, 5} };

  int face_neigh_elem[4][3] = {
                                {1, -1,-1},
                                {3, 0, 2},
                                {-1, 1, -1},
                                {-1, -1, 1} };
  /* *INDENT-ON* */

  if (cmesh == NULL) {
    /* If the cmesh is NULL. */
    return 0;
  }
  else {
    /* Checks if the cmesh was comitted. */
    retval = t8_cmesh_is_committed (cmesh);
    SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
    /* Checks for face consistency. */
    retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
    SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
    /* Checks if the number of elements was read correctly. */
    num_gtree = t8_cmesh_get_num_trees (cmesh);
    SC_CHECK_ABORT (num_gtree == number_elements,
                    "Number of elements in msh-file was read incorrectly.");
    /* Number of local trees. */
    lnum_trees = t8_cmesh_get_num_local_trees (cmesh);
    /* Iterate through the local elements and check if they were read properly. */
    /* All trees should be local to the master rank. */
    for (ltree_it = 0; ltree_it < lnum_trees; ltree_it++) {
      tree_class = t8_cmesh_get_tree_class (cmesh, ltree_it);
      check_element_type = t8_eclass_compare (tree_class, elem_type);
      SC_CHECK_ABORT (check_element_type == 0,
                      "Element type in msh-file was read incorrectly.");
      /* Get pointer to the vertices of the tree. */
      vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_it);
      /* Checking the msh-files elements and nodes. */
      for (int i = 0; i < 3; i++) {
        /* Checks if x and y coordinate of the nodes are not read correctly. */
        if (!((vertex[elements[ltree_it][i]][0] == (int) vertices[3 * i])
              && (vertex[elements[ltree_it][i]][1] ==
                  (int) vertices[(3 * i) + 1]))) {
          read_node = 0;
          SC_CHECK_ABORT (read_node == 1, "Node was read incorrectly.");
          /* Return error code, if the nodes are not read correctly. */
          return -1;
        }
        /* Checks whether the face neighbor elements are not read correctly. */
        ltree_id =
          t8_cmesh_get_face_neighbor (cmesh, ltree_it, i, NULL, NULL);
        if (!(ltree_id == face_neigh_elem[ltree_it][i])) {
          check_neigh_elem = 0;
          SC_CHECK_ABORT (check_neigh_elem == 1,
                          "The face neigbhor element in the example test file was not read correctly.");
          /* Return error code, if the face neighbor elements are not read correctly. */
          return -1;
        }
      }
    }
    /* If the checks were performed correctly. */
    return 1;
  }
}

/* Read the version 2 ascii file. This should work. */
static void
t8_test_cmesh_readmshfile_version2_ascii ()
{
  int                 retval;
  t8_cmesh_t          cmesh;
  const char          fileprefix[BUFSIZ - 4] =
    "test/testfiles/test_msh_file_vers2_ascii";
  char                filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_global_productionf ("Checking msh file version 2 ascii...\n");

  /* Check if file exists. */
  SC_CHECK_ABORTF (access (filename, R_OK) == 0, "Could not open file %s.\n",
                   filename);

  /* Try to read cmesh */
  cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0);
  SC_CHECK_ABORT (cmesh != NULL,
                  "Could not read cmesh from ascii version 2, but should be able to.");
  retval = t8_test_supported_msh_file (cmesh);
  SC_CHECK_ABORT (retval == 1, "Cmesh incorrectly read from file.");

  /* The cmesh was read sucessfully and we need to destroy it. */
  t8_cmesh_destroy (&cmesh);

  t8_global_productionf ("Could successfully read.\n");
}

/* Read version 2 binary file. We expect this to fail. */
static void
t8_test_cmesh_readmshfile_version2_bin ()
{
  int                 retval;
  t8_cmesh_t          cmesh;
  const char          fileprefix[BUFSIZ - 4] =
    "test/testfiles/test_msh_file_vers2_bin";
  char                filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_global_productionf ("Checking msh file version 2 binary...\n");

  /* Check if file exists. */
  SC_CHECK_ABORTF (access (filename, R_OK) == 0, "Could not open file %s.\n",
                   filename);

  /* Try to read cmesh */
  cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0);
  SC_CHECK_ABORT (cmesh == NULL,
                  "Expected fail of reading binary msh file v.2, but did not fail.");
  retval = t8_test_supported_msh_file (cmesh);
  SC_CHECK_ABORT (retval == 0,
                  "Unexpected return from t8_test_supported_msh_file.");

  t8_global_productionf ("Error handling successfull.\n");
}

/* Read version 4 ascii file. We expect this to fail. */
static void
t8_test_cmesh_readmshfile_version4_ascii ()
{
  int                 retval;
  t8_cmesh_t          cmesh;
  const char          fileprefix[BUFSIZ - 4] =
    "test/testfiles/test_msh_file_vers4_ascii";
  char                filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_global_productionf ("Checking msh file version 4 ascii...\n");

  /* Check if file exists. */
  SC_CHECK_ABORTF (access (filename, R_OK) == 0, "Could not open file %s.\n",
                   filename);

  /* Try to read cmesh */
  cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0);
  SC_CHECK_ABORT (cmesh == NULL,
                  "Expected fail of reading ascii msh file v4, but did not fail.");
  retval = t8_test_supported_msh_file (cmesh);
  SC_CHECK_ABORT (retval == 0,
                  "Unexpected return from t8_test_supported_msh_file.");
  t8_global_productionf ("Error handling successfull.\n");
}

/* Read version 4 bin file. We expect this to fail. */
static void
t8_test_cmesh_readmshfile_version4_bin ()
{
  int                 retval;
  t8_cmesh_t          cmesh;
  const char          fileprefix[BUFSIZ - 4] =
    "test/testfiles/test_msh_file_vers4_bin";
  char                filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix);

  t8_global_productionf ("Checking msh file version 4 binary...\n");

  /* Check if file exists. */
  SC_CHECK_ABORTF (access (filename, R_OK) == 0, "Could not open file %s.\n",
                   filename);

  /* Try to read cmesh */
  cmesh = t8_cmesh_from_msh_file (fileprefix, 1, sc_MPI_COMM_WORLD, 2, 0);
  SC_CHECK_ABORT (cmesh == NULL,
                  "Expected fail of reading binary msh file v4, but did not fail.");
  retval = t8_test_supported_msh_file (cmesh);
  SC_CHECK_ABORT (retval == 0,
                  "Unexpected return from t8_test_supported_msh_file.");
  t8_global_productionf ("Error handling successfull.\n");
}

int
main (int argc, char **argv)
{

  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

/* The mesh is unpartitioned and all trees are local to the master rank, no other rank has any trees. */

  t8_global_productionf ("Testing: Reading of msh-files.\n");

  /* Testing supported msh-file version 2 (as example). */
  t8_test_cmesh_readmshfile_version2_ascii ();

  /* Testing unsupported msh-file version 2 binary. */
  t8_test_cmesh_readmshfile_version2_bin ();

  /* Testing unsupported msh-file version 4 ascii. */
  t8_test_cmesh_readmshfile_version4_ascii ();

  /* Testing unsupported msh-file version 4 binary. */
  t8_test_cmesh_readmshfile_version4_bin ();

  t8_debugf ("Test successfull\n");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
