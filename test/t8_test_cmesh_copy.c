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

#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"


sc_MPI_Comm  comm_list[]={sc_MPI_COMM_WORLD,sc_MPI_COMM_NULL,sc_MPI_COMM_SELF};
size_t num_comm = sizeof(comm_list) / sizeof(comm_list[0]); 
int dim=3
int min_dim=1
int binary=2
int max_num_of_trees=100
int max_num_of_prisms=100

/* Test if a cmesh is committed properly and perform the
 * face consistency check. */
 
static void
test_cmesh_committed (t8_cmesh_t cmesh)
{
  int                 retval;

  retval = t8_cmesh_is_committed (cmesh);
  SC_CHECK_ABORT (retval == 1, "Cmesh commit failed.");
  retval = t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees);
  SC_CHECK_ABORT (retval == 1, "Cmesh face consistency failed.");
}


/*
 *The functions t8_get_*_cmesh_testcases return the number of 
 *testcases for a given cmesh type.
 */
static  int  t8_get_comm_only_cmesh_testcases():
{
  /*number of functions that only take comm as input * number of comm */
  return 10*num_comm;
}

static  int  t8_get_new_hypercube_cmesh_testcases():
{
  /* number of element types * number of comm*3 binary variables(do_bcast, do_partition,periodic) */
  return T8_ECLASS_COUNT*num_comm*binary*binary*binary;
}

static  int  t8_get_new_empty_cmesh_testcases():
{
  /* number of comm* 1 binary variable(do_partition)* possible dimensions(check dim 0 to 4, therefore we need to add 2 to dim)*/
  return num_comm*binary*(dim+2);
}

static  int  t8_get_new_from_class_cmesh_testcases():
{
  /* number of element types * number of comm */
  return T8_ECLASS_COUNT*num_comm;
}

static  int  t8_get_new_hypercube_hybrid_cmesh_testcases():
{
  /* possible dim * number of comm*2 binary variables(do_partition,periodic) */
  return 4*num_comm*binary*binary;
}

static  int t8_get_new_periodic_cmesh_testcases():
{
  /*  number of comm * possible dim */
  return num_comm*dim;
}

static  int  t8_get_new_bigmesh_cmesh_testcases():
{
  /*  number of element types * number of trees * number of comm */
  return T8_ECLASS_COUNT*max_num_of_trees*num_comm;
}

static  int  t8_get_new_prism_cake_cmesh_testcases():
{
  /*  number of comm * number of prisms */
  return num_comm*max_num_of_prisms;
}

static  int t8_get_new_disjoint_bricks_cmesh_testcases():
{
  /*  TO DO: add comment*/
  return 20*20*20*binary*binary*binary*num_comm;
}

static int t8_get_all_testcases():
{
  /* The number of all tests.*/
  return t8_get_comm_only_cmesh_testcases()+t8_get_new_hypercube_cmesh_testcases()+t8_get_new_empty_cmesh_testcases()+t8_get_new_from_class_cmesh_testcases()+t8_get_new_hypercube_hybrid_cmesh_testcases()+t8_get_new_periodic_cmesh_testcases()+t8_get_new_bigmesh_cmesh_testcases()+t8_get_new_prism_cake_cmesh_testcases()+t8_get_new_disjoint_bricks_cmesh_testcases();
}


static              t8_cmesh_t
t8_test_create_comm_only_cmesh(int cmesh_id)
{
  switch ((int)(cmesh_id/3)) {
  case 0:
    return t8_cmesh_new_periodic_tri (comm[cmesh_id % 3]);
  case 1:
    return t8_cmesh_new_periodic_hybrid (comm[cmesh_id % 3]);
  case 2:
    return t8_cmesh_new_periodic_line_more_trees (comm[cmesh_id % 3]);
  case 3:
    return t8_cmesh_new_line_zigzag (comm[cmesh_id % 3]);
  case 4:
    return t8_cmesh_new_prism_deformed (comm[cmesh_id % 3]);
  case 5:
    return t8_cmesh_new_prism_cake_funny_oriented (comm[cmesh_id % 3]);
  case 6:
    return t8_cmesh_new_prism_geometry (comm[cmesh_id % 3]);
  case 7:
    return t8_cmesh_new_tet_orientation_test (comm[cmesh_id % 3]);
  case 8:
    return t8_cmesh_new_hybrid_gate (comm[cmesh_id % 3]);
  case 9:
    return t8_cmesh_new_hybrid_gate_deformed (comm[cmesh_id % 3]);
  }
}

static              t8_cmesh_t
t8_test_create_new_hypercube_cmesh(int cmesh_id)
{
  return t8_cmesh_new_hypercube ((t8_eclass_t) (int)(cmesh_id/(num_comm*binary*binary*binary) %(num_comm*binary*binary*binary)),
                                 comm[((int)cmesh_id/binary*binary*binary)%binary*binary], (int)((cmesh_id%binary*binary*binary)/binary*binary),
                                 ((int)(cmesh_id/binary))%binary, cmesh_id % binary);
}

static              t8_cmesh_t
t8_test_create_new_empty_cmesh(int cmesh_id)
{
  return t8_cmesh_new_empty ((int)((cmesh_id/8)%8), (int)((cmesh_id/(dim+2))%binary), cmesh_id % (dim+2));
}

static              t8_cmesh_t
t8_test_create_new_from_class_cmesh(int cmesh_id)
{
  return t8_cmesh_new_from_class ((t8_eclass_t)((int) (cmesh_id/3 % T8_ECLASS_COUNT)), comm[cmesh_id % 3]);
}

static              t8_cmesh_t
t8_test_create_new_hypercube_hybrid_cmesh(int cmesh_id)
{
  return t8_cmesh_new_hypercube_hybrid (((int)((cmesh_id/12)%12))%4, comm[(int)((cmesh_id/3)%3)], ((int)(cmesh_id/2))%2, cmesh_id%2);
}


static              t8_cmesh_t
t8_test_create_new_periodic_cmesh(int cmesh_id)
{
  return t8_cmesh_new_periodic (comm[(int)((cmesh_id/4)%8)], (cmesh_id%dim)+min_dim);
}

static              t8_cmesh_t
t8_test_create_new_bigmesh_cmesh(int cmesh_id)
{
  return t8_cmesh_new_bigmesh ((t8_eclass_t)((int)(cmesh_id/300)%300), (int) ((cmesh_id/3)%100), comm[cmesh_id%2]);
}

static              t8_cmesh_t
t8_test_create_new_prism_cake_cmesh(int cmesh_id)
{
  return t8_cmesh_new_prism_cake (comm[(int)((cmesh_id/100)%3)], (cmesh_id% 100)+3);
}

static              t8_cmesh_t
t8_test_create_new_disjoint_bricks_cmesh(int cmesh_id)
{
  return t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x, t8_gloidx_t num_y, t8_gloidx_t num_z,int x_periodic,int y_periodic,int z_periodic,sc_MPI_Comm comm);
}

static              t8_cmesh_t
t8_test_create_cmesh(int cmesh_id)
{
  if(0<=cmesh_id<t8_get_comm_only_cmesh_testcases()){
    return t8_test_create_comm_only_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()<t8_get_new_hypercube_cmesh_testcases()){
    return t8_test_create_new_hypercube_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases() 
                      - t8_get_new_hypercube_cmesh_testcases()<t8_get_new_empty_cmesh_testcases()){
    return t8_test_create_new_empty_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()<t8_get_new_from_class_cmesh_testcases()){
    return t8_test_create_new_from_class_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()
                      - t8_get_new_from_class_cmesh_testcases()<t8_get_new_hypercube_hybrid_cmesh_testcases()){
    return t8_test_create_new_hypercube_hybrid_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()
                      - t8_get_new_from_class_cmesh_testcases()
                      - t8_get_new_hypercube_hybrid_cmesh_testcases()<t8_get_new_periodic_cmesh_testcases()){
    return t8_test_create_new_periodic_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()
                      - t8_get_new_from_class_cmesh_testcases()
                      - t8_get_new_hypercube_hybrid_cmesh_testcases()
                      - t8_get_new_periodic_cmesh_testcases()<t8_get_new_bigmesh_cmesh_testcases()){
    return t8_test_create_new_bigmesh_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()
                      - t8_get_new_from_class_cmesh_testcases()
                      - t8_get_new_hypercube_hybrid_cmesh_testcases()
                      - t8_get_new_periodic_cmesh_testcases()
                      - t8_get_new_bigmesh_cmesh_testcases()<t8_get_new_prism_cake_cmesh_testcases()){
    return t8_test_create_new_prism_cake_cmesh(cmesh_id);
  }
  else if(0<=cmesh_id - t8_get_comm_only_cmesh_testcases()
                      - t8_get_new_hypercube_cmesh_testcases()
                      - t8_get_new_empty_cmesh_testcases()
                      - t8_get_new_from_class_cmesh_testcases()
                      - t8_get_new_hypercube_hybrid_cmesh_testcases()
                      - t8_get_new_periodic_cmesh_testcases()
                      - t8_get_new_bigmesh_cmesh_testcases()
                      - t8_get_new_prism_cake_cmesh_testcases()<t8_get_new_disjoint_bricks_cmesh_testcases()){
    return t8_test_create_new_disjoint_bricks_cmesh(cmesh_id);
  }
}

static void
test_cmesh_copy (int cmesh_id,comm)
{
  int                 retval;
  t8_cmesh_t          cmesh_original, cmesh_copy;

      /* Create new cmesh */
      cmesh_original = t8_test_create_cmesh (cmesh_id);
      test_cmesh_committed (cmesh_original);
      /* Set up the cmesh copy */
      t8_cmesh_init (&cmesh_copy);
      /* We need the original cmesh later, so we ref it */
      t8_cmesh_ref (cmesh_original);
      t8_cmesh_set_derive (cmesh_copy, cmesh_original);
      /* Commit and check commit */
      t8_cmesh_commit (cmesh_copy, comm);
      test_cmesh_committed (cmesh_copy);
      /* Check for equality */
      retval = t8_cmesh_is_equal (cmesh_copy, cmesh_original);
      SC_CHECK_ABORT (retval == 1, "Cmesh copy failed.");
      /* Clean-up */
      t8_cmesh_destroy (&cmesh_copy);
      t8_cmesh_destroy (&cmesh_original);
}

static void
test_cmesh_copy_all(comm)
{
  /* Test all cmeshes over all different inputs we get through their id */
  for (cmesh_id = 0; cmesh_id < t8_get_all_testcases(); cmesh_id++) {
    test_cmesh_copy (cmesh_id,comm);
  }
}


int
main (int argc, char **argv)
{
  int                 mpiret, i, eci;
  sc_MPI_Comm         comm;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  comm = sc_MPI_COMM_WORLD;
  sc_init (comm, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf ("Testing cmesh copy.\n");
  test_cmesh_copy_all(comm);
  t8_global_productionf ("Done testing cmesh copy.\n");
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
