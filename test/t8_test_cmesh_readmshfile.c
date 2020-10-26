//#include <t8_eclass.h>
#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include <t8_cmesh_readmshfile.h>

/* Define constants for the properties of the example test case of the tested files. */
#define T8_TEST_READING_MSH_FILE_DIM 2
#define T8_TEST_READING_MSH_FILE_NUM_ELEM 3
#define T8_TEST_READING_MSH_FILE_ELEM_TYPE 3

#if 1
static int
t8_test_supported_msh_file (t8_cmesh_t cmesh)
{

  int                 retval;
  int                 read_node;
  double             *vertices;
  t8_gloidx_t         num_gtree;
  //t8_locidx_t ltree_id;
  t8_locidx_t         ltree_it;
  t8_locidx_t         lnum_trees;
  t8_eclass_t         tree_class;
  //t8_ctree_t ltree;
  //double vertex[6][2] = {{0.5,1},{1,0},{1.5,1},{0,0},{2,0},{1,2}};
  //int elements[4][3] = {{0,1,2},{0,3,1},{1,4,2},{0,2,5}};
  double              vertex[6][2] =
    { {0, 0}, {0, 1}, {0, 2}, {0.5, 1}, {1.5, 1}, {1, 2} };
  int                 elements[4][3] =
    { {0, 1, 3}, {1, 4, 3}, {1, 2, 4}, {3, 4, 5} };

  if (cmesh == NULL) {
    return 0;
  }
  else {
    /* Checks if the cmesh was comitted. */
    retval = t8_cmesh_is committed (cmesh);
    SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
    /* Checks for face consistency. */
    retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
    SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
    /* Checks if the number of elements was read correctly. */
    num_gtree = t8_cmesh_get_num_trees (cmesh);
    SC_CHECK_ABORT (num_gtree == T8_TEST_READING_MSH_FILE_NUM_ELEM,
                    "Number of elements in msh-file was read incorrectly.");
    /* Number of local trees. */
    lnum_trees = t8_cmesh_get_num_local_trees (cmesh);
    /* Iterate through the local elements and check if they were read properly. */
    /* All trees should be local to the master rank. */
    for (ltree_it = 0; ltree_it < lnum_trees; ltree_it++) {
      tree_class = t8_cmesh_get_tree_class (cmesh, ltree_it);
      SC_CHECK_ABORT (tree_class == T8_TEST_READING_MSH_FILE_ELEM_TYPE,
                      "Element type in msh-file was read incorrectly.");
      /* Get pointer to the vertices of the tree. */
      *vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_it);
      /* Checking the nodes of the element. */
      read_node = 0;
      for (int i = 0; i < 3, i++) {
        /* Checks if x and y coordinate are read correctly. */
        if (vertex[elements[ltree_it][i]][0] == vertices[i * 3]
            && vertex[elements[ltree_it][i]][1] == vertices[i * 3 + 1]) {
          read_node = 1;
        }
        else {
          read_node = 0;
        }
      }
      /* Case if the node does not belong the element, and therefore, must have been read wrong. */
      SC_CHECK_ABORT (read_node == 1, "Node was read incorrectly.");
    }
    /* Check for face neighbors. */
  }
}

#endif

int
main (int argc, char **argv)
{

  int                 mpiret;
  int                 partition, dim, master, retval;
  t8_cmesh_t          cmesh;

/* Declaration of test variables. */
  partition = 1;
  dim = 2;
  master = 0;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

/* The mesh is unpartitioned and all trees are local to the master rank, no other rank has any trees. */

  t8_debugf ("Testing: Reading of msh-files.\n");
#if 0
  /* Testing supported msh-file version (as example). */

  cmesh =
    t8_cmesh_from_msh_file ("./test_msh_file_vers2_ascii", partition,
                            sc_MPI_COMM_WORLD, dim, master);
  retval = t8_test_supported_msh_file (cmesh);
  if (!retval) {
    t8_global_errorf ("Reading of supported msh-file version failed.\n");
    goto end_test;
  }
  else {
    t8_global_productionf
      ("Test: Reading of msh-files of supported version was succesful.\n");
  }
#endif
  t8_debugf ("....Test MSH-File wurde aufgerufen....\n");
#if 0
  /* Testing unsupported version of msh-files (bin-format). */
  cmesh =
    t8_cmesh_from_msh_file ("./test_msh_file_vers2_bin", partition,
                            sc_MPI_COMM_WORLD, dim, master);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval) {
    t8_global_errorf ("Rejecting of unsupported msh-files failed.\n");
    goto end_test;
  }

  /* Testing unsupported version of msh-files (version 4). */
  cmesh =
    t8_cmesh_from_msh_file ("./test_msh_file_vers4_ascii", partition,
                            sc_MPI_COMM_WORLD, dim, master);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval) {
    t8_global_errorf ("Rejecting of unsupported msh-files failed.\n");
    goto end_test;
  }

  /* Testing unsupported version of msh-files (version 1). */
  cmesh =
    t8_cmesh_from_msh_file ("./test_msh_file_vers1_ascii", partition,
                            sc_MPI_COMM_WORLD, dim, master);
  retval = t8_test_supported_msh_file (cmesh);
  if (retval) {
    t8_global_errorf ("Rejecting of unsupported msh-files failed.\n");
    goto end_test;
  }
#endif

  t8_cmesh_destroy (&cmesh);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;

/* In case anything failed, the test is going to be closed. */
end_test:
  t8_cmesh_destroy (&cmesh);
  return 0;
}
