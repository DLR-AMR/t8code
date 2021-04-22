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

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_partition.h"
#include <t8_eclass.h>

#define T8_CMESH_TEST_NUM_COMMS 2
#define T8_CMESH_BINARY 2
#define T8_CMESH_DIM_RANGE 4    /* this is the dim range for hypercube hybrid and empty cmesh */
#define T8_CMESH_MAX_TEST_DIMS 3
#define T8_CMESH_MAX_NUM_OF_TREES 100
#define T8_CMESH_MAX_NUM_OF_PRISMS 100
#define T8_CMESH_MAX_NUM_XYZ_TREES 20

/* The functions t8_get_*_cmesh_testcases return the number of 
 * testcases for a given cmesh type. To get all possible inputs, we multiply 
 * the number of different inputs of all variables. The function t8_get_comm_only_cmesh_testcases()
 * returns the number of testcases for all cmeshes that only take a comm as input all added together. */
int
t8_get_number_of_comm_only_cmesh_testcases ()
{
  const int           num_of_comm_only_func = 10;
  /* Number of testcases = Number of functions that only take communicator as input * Number of communicators */
  return num_of_comm_only_func * T8_CMESH_TEST_NUM_COMMS;
}

/* The function t8_get_new_hypercube_cmesh_testcases() returns the number of testcases for 
 * the hypercube cmesh. */
int
t8_get_number_of_new_hypercube_cmesh_testcases ()
{
  /* Number of testcases = number of element types * number of comm * 
   * 3 T8_CMESH_BINARY variables(do_bcast, do_partition,periodic) */
  return T8_ECLASS_COUNT * T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY *
    T8_CMESH_BINARY * T8_CMESH_BINARY;
}

/* The function t8_get_new_empty_cmesh_testcases() returns the number of testcases for 
 * the empty cmesh. */
int
t8_get_number_of_new_empty_cmesh_testcases ()
{
  /* Number of testcases = number of comm * 1 binary (do_partition)* possible dimensions(check dim 0 to 4) */
  return T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY * T8_CMESH_DIM_RANGE;
}

/* The function t8_get_new_from_class_cmesh_testcases() returns the number of testcases for 
 * the new_from_class cmesh. */
int
t8_get_number_of_new_from_class_cmesh_testcases ()
{
  /* Number of testcases = number of element types * number of comm */
  return T8_ECLASS_COUNT * T8_CMESH_TEST_NUM_COMMS;
}

/* The function t8_get_new_hypercube_hybrid_cmesh_testcases() returns the number of testcases for 
 * the new_hypercube_hybrid cmesh. */
int
t8_get_number_of_new_hypercube_hybrid_cmesh_testcases ()
{
  /* Number of testcases = possible dim * number of comm*2 binary variables(do_partition,periodic) */
  return T8_CMESH_MAX_TEST_DIMS * T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY *
    T8_CMESH_BINARY;
}

/* The function t8_get_new_periodic_cmesh_testcases() returns the number of testcases for 
 * the new_periodic cmesh. */
int
t8_get_number_of_new_periodic_cmesh_testcases ()
{
  /* Number of testcases = number of comm * possible dim */
  return T8_CMESH_TEST_NUM_COMMS * T8_CMESH_MAX_TEST_DIMS;
}

/* The function t8_get_new_bigmesh_cmesh_testcases() returns the number of testcases for 
 * the new_bigmesh cmesh. */
int
t8_get_number_of_new_bigmesh_cmesh_testcases ()
{
  /* Number of testcases = number of element types * number of trees * number of comm */
  return T8_ECLASS_COUNT * T8_CMESH_MAX_NUM_OF_TREES *
    T8_CMESH_TEST_NUM_COMMS;
}

/* The function t8_get_new_prism_cake_cmesh_testcases() returns the number of testcases for 
 * the new_prism_cake cmesh. */
int
t8_get_number_of_new_prism_cake_cmesh_testcases ()
{
  /* Number of testcases = number of comm * number of prisms */
  return T8_CMESH_TEST_NUM_COMMS * T8_CMESH_MAX_NUM_OF_PRISMS;
}

/* The function t8_get_new_disjoint_bricks_cmesh_testcases() returns the number of testcases for 
 * the new_disjoint_bricks cmesh. */
int
t8_get_number_of_new_disjoint_bricks_cmesh_testcases ()
{
  /* Number of testcases = num trees in x direction * num trees in y direction * num trees in z direction 
   * 3 binary variables (x_periodic, y_periodic, z_periodic)*num of comm */
  return T8_CMESH_MAX_NUM_XYZ_TREES * T8_CMESH_MAX_NUM_XYZ_TREES *
    T8_CMESH_MAX_NUM_XYZ_TREES * T8_CMESH_BINARY * T8_CMESH_BINARY *
    T8_CMESH_BINARY * T8_CMESH_TEST_NUM_COMMS;
}

/* The function t8_get_all_testcases() returns the number of testcases for 
 * all cmesh types. We need to know this, because test_cmesh_copy_all needs to know how 
 * many ids to check. */
int
t8_get_number_of_all_testcases ()
{
  /* The number of all tests for all cmesh types is the sum of all individual testcase numbers. */
  return (t8_get_number_of_comm_only_cmesh_testcases () +
          t8_get_number_of_new_hypercube_cmesh_testcases ()
          + t8_get_number_of_new_empty_cmesh_testcases () +
          t8_get_number_of_new_from_class_cmesh_testcases ()
          + t8_get_number_of_new_hypercube_hybrid_cmesh_testcases () +
          t8_get_number_of_new_periodic_cmesh_testcases ()
          + t8_get_number_of_new_bigmesh_cmesh_testcases () +
          t8_get_number_of_new_prism_cake_cmesh_testcases ()
          + t8_get_number_of_new_disjoint_bricks_cmesh_testcases ());
}
