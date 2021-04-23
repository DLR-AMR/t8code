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
#include <t8_eclass.h>
#include "t8_cmesh/t8_cmesh_testcases.h"

#define T8_CMESH_TEST_NUM_COMMS 2
#define T8_CMESH_BINARY 2
#define T8_CMESH_DIM_RANGE 4    /* this is the dim range for hypercube hybrid and empty cmesh */
#define T8_CMESH_MAX_TEST_DIMS 3
#define T8_CMESH_MAX_NUM_OF_TREES 100
#define T8_CMESH_MAX_NUM_OF_PRISMS 100
#define T8_CMESH_MAX_NUM_XYZ_TREES 20

sc_MPI_Comm         t8_comm_list[T8_CMESH_TEST_NUM_COMMS] =
  { sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF };
char                t8_comm_string_list[2][18] =
  { "sc_MPI_COMM_WORLD", "sc_MPI_COMM_SELF" };

/* The functions t8_test_create_*_cmesh create a cmesh of a given type with a unique input depending on the cmesh_id. 
 * Every cmesh-type has a certain number of input variables. We're interested in testing cmeshes with different inputs.
 * For this we need a way of running through all different inputs. In t8_cmesh_testcases.c we have functions that return
 * the number of inputs we have to check for every cmesh. In this file we actually create the cmeshes. We use a consistent
 * way of running over the variables, e.g. t8_cmesh_new_hypercube (eclass, comm, do_bcast, do_partition, periodic):
 * eclass(0 to T8_ECLASS_COUNT-1) |comm(0 to T8_CMESH_TEST_NUM_COMMS-1)| do_bcast | do_partition | periodic (the last 3 are binary)
 *              0                 |                0                   |    0     |       0      |             0
 *              0                 |                0                   |    0     |       0      |             1
 *              0                 |                0                   |    0     |       1      |             0
 *              0                 |                0                   |    0     |       1      |             1
 *              0                 |                0                   |    1     |       0      |             0
 *              0                 |                0                   |    1     |       0      |             1
 *              0                 |                0                   |    1     |       1      |             0
 *              0                 |                0                   |    1     |       1      |             1
 *              0                 |                1                   |    0     |       0      |             0
 *             ...                |               ...                  |   ...    |      ...     |            ...
 *              7                 |                1                   |    1     |       1      |             1
 * The idea is to get these sequences by using a cmesh_id, meaning e.g. cmesh_id=0 gives us the first row,
 * cmesh_id=1 the second row and so on. In each function you will find lines like:
 * do_partition = (cmesh_id / T8_CMESH_BINARY) % T8_CMESH_BINARY (*)
 * We have to understand that the division is integer division and therefore like taking the floor of the division. 
 * The modulo makes sure that we are in the right range. A binary variable should only be between 0 and 1. 
 * To understand equation (*), consider the following: Since we use integer division by 2, the bracket gives us 
 * 0 , 1/2 , 1 , 3/2 , 2 , 5/2... taking the floor gives: 0,0,1,1,2,2,3,3,...
 * Taking modulo 2 turns this into: 0,0,1,1,0,0,1,1,...
 * So the division is important so that we only change every x steps if we divide by x.
 * The modulo is so that we stay in the correct range. 
 */

/* The function t8_test_create_comm_only_cmesh(int cmesh_id) returns the wanted cmesh with the wanted comm for the given id. 
 * The comm is taken from the t8_comm_list. The switch inside t8_test_create_comm_only_cmesh(int cmesh_id)
 * chooses the cmesh-type. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique cmesh which only takes comm as input.
 */
t8_cmesh_t
t8_test_create_comm_only_cmesh (int cmesh_id)
{
  /* Because we get i  through integer division, it is essentially the floor.
   * The variable i gives us the cmesh type because we use the switch. It switches
   * every T8_CMESH_TEST_NUM_COMMS steps(e.g. for T8_CMESH_TEST_NUM_COMMS=2 0,0,1,1,2,2,..). 
   * The comm changes every step and is in range of the number of comm we have. 
   * The modulo does exactly this. The comm_num changes with increasing cmesh_id as (0,1,0,1,0,1,0,1,...). 
   * This makes it possible to check every cmesh that takes only comm with every comm in t8_comm_list. 
   * Notice that each combination is unique.
   */
  int                 i = cmesh_id / T8_CMESH_TEST_NUM_COMMS;
  int                 comm_num = cmesh_id % T8_CMESH_TEST_NUM_COMMS;
  switch (i) {
  case 0:
    t8_debugf ("Creating new periodic tri cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_periodic_tri (t8_comm_list[comm_num]);
  case 1:
    t8_debugf ("Creating new periodic hybrid cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_periodic_hybrid (t8_comm_list[comm_num]);
  case 2:
    t8_debugf ("Creating new periodic line more trees cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_periodic_line_more_trees (t8_comm_list[comm_num]);
  case 3:
    t8_debugf ("Creating new line zig zag cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_line_zigzag (t8_comm_list[comm_num]);
  case 4:
    t8_debugf ("Creating new prism deformed cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_prism_deformed (t8_comm_list[comm_num]);
  case 5:
    t8_debugf ("Creating new prism cake funny oriented cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_prism_cake_funny_oriented (t8_comm_list[comm_num]);
  case 6:
    t8_debugf ("Creating new prism geometry cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_prism_geometry (t8_comm_list[comm_num]);
  case 7:
    t8_debugf ("Creating new tet orientation test cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_tet_orientation_test (t8_comm_list[comm_num]);
  case 8:
    t8_debugf ("Creating new hybrid gate cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_hybrid_gate (t8_comm_list[comm_num]);
  case 9:
    t8_debugf ("Creating new hybrid gate deformed cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_hybrid_gate_deformed (t8_comm_list[comm_num]);
  default:
    t8_debugf ("Creating new periodic tri cmesh. comm=%s \n",
               t8_comm_string_list[comm_num]);
    return t8_cmesh_new_periodic_tri (t8_comm_list[comm_num]);
  }
}

/* The function t8_test_create_new_hypercube_cmesh(int cmesh_id) returns a new hypercube cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. It avoids the case (eclass = pyramid & periodic=1) since this is not allowed. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_hypercube cmesh.
 */
t8_cmesh_t
t8_test_create_new_hypercube_cmesh (int cmesh_id)
{
  /* The variable periodic changes every step between 0 and 1(so 0,1,0,1,...), because these are all 
   * cases for a binary variable. do_partition is also binary, but should change whenever periodic 
   * looped one time through its values, so it is changing only every second step(0,0,1,1,0,0,1,1,...).
   * do_bcast should as a binary variable then change every 4 steps(0,0,0,0,1,1,1,1,0,...),
   * because then do_partition looped once. eclass should be between 0 and T8_ECLASS_COUNT, 
   * therefore we take modulo T8_ECLASS_COUNT. It should change every 16 steps because then the 
   * comm looped once. This way for every cmesh_id we get a unique combination of input variables.
   */
  int                 comm_num = (cmesh_id
                                  / (T8_CMESH_BINARY * T8_CMESH_BINARY *
                                     T8_CMESH_BINARY))
    % (T8_CMESH_TEST_NUM_COMMS);
  t8_eclass_t         eclass = (cmesh_id
                                / (T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY *
                                   T8_CMESH_BINARY * T8_CMESH_BINARY))
    % T8_ECLASS_COUNT;
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 do_bcast = (cmesh_id
                                  % (T8_CMESH_BINARY * T8_CMESH_BINARY *
                                     T8_CMESH_BINARY))
    / (T8_CMESH_BINARY * T8_CMESH_BINARY);
  int                 do_partition =
    (cmesh_id / T8_CMESH_BINARY) % T8_CMESH_BINARY;
  int                 periodic = cmesh_id % T8_CMESH_BINARY;

  t8_debugf
    ("Creating new hypercube cmesh. eclass=%s, comm=%s, do_bcast=%i, do_partition=%i, periodic=%i \n",
     t8_eclass_to_string[eclass], t8_comm_string_list[comm_num], do_bcast,
     do_partition, periodic);

  if ((int) eclass == 7 && periodic == 1) {
    return t8_cmesh_new_hypercube (eclass, comm, do_bcast, do_partition, 0);
  }

  return t8_cmesh_new_hypercube (eclass, comm, do_bcast, do_partition,
                                 periodic);
}

/* The function t8_test_create_new_empty_cmesh(int cmesh_id) returns a new empty cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_empty cmesh.
 */
t8_cmesh_t
t8_test_create_new_empty_cmesh (int cmesh_id)
{
  /* The variable dim changes every step between 0 and T8_CMESH_DIM_RANGE(so 0,1,2,3,0,1,...), because these are all 
   * cases for a the dim variable. do_partition is binary, but should change whenever dim 
   * looped one time through its values, so it is changing only every fourth step(0,0,0,0,1,1,1,1,0,..).
   * The comm should then change every 8 steps, because then do_partition loops once. It also needs to be 
   * in range of T8_CMESH_TEST_NUM_COMMS, that is why we take modulo T8_CMESH_TEST_NUM_COMMS.
   * This way we get a unique input for every cmesh_id.
   *            comm                |           do_partition             |   dim    
   *              0                 |                0                   |    0     
   *              0                 |                0                   |    1     
   *              0                 |                0                   |    2     
   *              0                 |                0                   |    3     
   *              0                 |                1                   |    0     
   *              0                 |                1                   |    1     
   *              0                 |                1                   |    2     
   *              0                 |                1                   |    3     
   *              1                 |                0                   |    0     
   *              1                 |                0                   |    1     
   *              1                 |                0                   |    2     
   *              1                 |                0                   |    3     
   *              1                 |                1                   |    0     
   *              1                 |                1                   |    1     
   *              1                 |                1                   |    2     
   *              1                 |                1                   |    3   
   */
  int                 comm_num = ((cmesh_id
                                   / (T8_CMESH_BINARY * T8_CMESH_DIM_RANGE))
                                  % T8_CMESH_TEST_NUM_COMMS);
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 do_partition =
    (cmesh_id / T8_CMESH_DIM_RANGE) % T8_CMESH_BINARY;
  int                 dim = cmesh_id % T8_CMESH_DIM_RANGE;

  t8_debugf ("Creating new empty cmesh. comm=%s , do_partition=%i, dim=%i \n",
             t8_comm_string_list[comm_num], do_partition, dim);

  return t8_cmesh_new_empty (comm, do_partition, dim);
}

/* The function t8_test_create_new_from_class_cmesh(int cmesh_id) returns a new create_new_from_class cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_from_class cmesh.
 */
t8_cmesh_t
t8_test_create_new_from_class_cmesh (int cmesh_id)
{
 /* The variable comm_num changes every step between 0 and T8_CMESH_TEST_NUM_COMMS(so 0,1,0,1,0,1,...), because these are all 
  * cases for the comm_num variable. eclass should then change every T8_CMESH_TEST_NUM_COMMS steps, because 
  * then comm_num looped once through all its values. That is why we dividie by T8_CMESH_TEST_NUM_COMMS.
  * It also must be in range of T8_ECLASS_COUNT, so we take modulo T8_ECLASS_COUNT. 
  * This way we get a unique input for every cmesh_id.
  *            eclass              |              comm             
  *              0                 |                0     
  *              0                 |                1     
  *              1                 |                0     
  *              1                 |                1     
  *              2                 |                0        
  *              2                 |                1       
  *              3                 |                0          
  *              3                 |                1       
  *              4                 |                0       
  *              4                 |                1        
  *              5                 |                0        
  *              5                 |                1         
  *              6                 |                0        
  *              6                 |                1         
  *              7                 |                0         
  *              7                 |                1         
  */
  int                 comm_num = cmesh_id % T8_CMESH_TEST_NUM_COMMS;
  t8_eclass_t         eclass =
    (cmesh_id / T8_CMESH_TEST_NUM_COMMS) % T8_ECLASS_COUNT;
  sc_MPI_Comm         comm = t8_comm_list[comm_num];

  t8_debugf ("Creating new_from_class cmesh. eclass=%s, comm=%s \n",
             t8_eclass_to_string[eclass], t8_comm_string_list[comm_num]);

  return t8_cmesh_new_from_class (eclass, comm);
}

/* The function t8_test_create_new_hypercube_hybrid_cmesh(int cmesh_id) returns a new_hypercube_hybrid cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_hypercube_hybrid cmesh. 
 */
t8_cmesh_t
t8_test_create_new_hypercube_hybrid_cmesh (int cmesh_id)
{
 /* The variable periodic changes every step between 0 and 1(so 0,1,0,1,0,1,...), because these are all 
  * cases for a binary variable. do_partition should then change every second step, because 
  * then periodic looped once through all its values. That is why we dividie by T8__CMESH_BINARY.
  * It also must be in range of T8__CMESH_BINARY, so we take modulo T8__CMESH_BINARY(0,0,1,1,0,...). 
  * comm_num should change every 4 steps, so we divide by 4 and for it to be in range 
  * of T8_CMESH_TEST_NUM_COMMS, we take modulo T8_CMESH_TEST_NUM_COMMS. dim should change 
  * every 8 steps, because then comm_num looped once through all its values,so we divide by 4 
  * For it to be in range of T8_CMESH_DIM_RANGE, we take modulo T8_CMESH_DIM_RANGE.
  * This way we get a unique input for every cmesh_id.
  *    dim     |     comm    | do_partition | periodic
  *     0      |       0     |       0      |     0
  *     0      |       0     |       0      |     1
  *     0      |       0     |       1      |     0
  *     0      |       0     |       1      |     1
  *     0      |       1     |       0      |     0
  *     0      |       1     |       0      |     1
  *     0      |       1     |       1      |     0
  *     0      |       1     |       1      |     1
  *     1      |       0     |       0      |     0
  *    ...     |      ...    |      ...     |    ...     
  *     3      |       1     |       1      |     1     
  */
  int                 comm_num = (cmesh_id
                                  /( T8_CMESH_BINARY * T8_CMESH_BINARY))
    % T8_CMESH_TEST_NUM_COMMS;
  int                 dim = (cmesh_id
                              / (T8_CMESH_BINARY * T8_CMESH_BINARY *
                                 T8_CMESH_TEST_NUM_COMMS))
                             % (T8_CMESH_BINARY * T8_CMESH_BINARY *
                                T8_CMESH_TEST_NUM_COMMS)
    % (T8_CMESH_DIM_RANGE);
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 do_partition =
    (cmesh_id / T8_CMESH_BINARY) % T8_CMESH_BINARY;
  int                 periodic = cmesh_id % T8_CMESH_BINARY;

  t8_debugf
    ("Creating new hypercube hybrid cmesh. dim=%i, comm=%s , do_partition=%i, periodic=%i \n",
     dim, t8_comm_string_list[comm_num], do_partition, periodic);

  return t8_cmesh_new_hypercube_hybrid (dim, comm, do_partition, periodic);
}

/* The function t8_test_create_new_periodic_cmesh(int cmesh_id) returns a new_periodic cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal dimension is 1.
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_periodic cmesh. 
 */
t8_cmesh_t
t8_test_create_new_periodic_cmesh (int cmesh_id)
{ 
 /*
  *    comm    |      dim    
  *     0      |       1     
  *     0      |       2     
  *     0      |       3    
  *     1      |       1     
  *     1      |       2    
  *     1      |       3    
  */
  const int           min_dim = 1;
  int                 comm_num = (cmesh_id / T8_CMESH_TEST_NUM_COMMS)
    % T8_CMESH_TEST_NUM_COMMS;
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 dim = (cmesh_id % T8_CMESH_MAX_TEST_DIMS) + min_dim;

  t8_debugf ("Creating new periodic cmesh. comm=%s,dim=%i \n",
             t8_comm_string_list[comm_num], dim);

  return t8_cmesh_new_periodic (comm, dim);
}

/* The function t8_test_create_new_bigmesh_cmesh(int cmesh_id) returns a new_bigmesh cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal number of trees is 1. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_bigmesh cmesh.
 */
t8_cmesh_t
t8_test_create_new_bigmesh_cmesh (int cmesh_id)
{
 /*
  *    eclass  |   num_trees |     comm   
  *     0      |       0     |       0 
  *     0      |       0     |       1  
  *     0      |       1     |       0  
  *     0      |       1     |       1  
  *     0      |       2     |       0  
  *     0      |       2     |       1    
  *     0      |       3     |       0     
  *     0      |       3     |       1    
  *     0      |       4     |       0    
  *    ...     |      ...    |      ...   
  *     7      |      100    |       1    
  */
  const int           min_num_trees = 1;
  int                 comm_num = cmesh_id % T8_CMESH_BINARY;
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 num_trees = min_num_trees
    + ((cmesh_id / T8_CMESH_TEST_NUM_COMMS)
       % T8_CMESH_MAX_NUM_OF_TREES);
  t8_eclass_t         eclass = (cmesh_id
                                / (T8_CMESH_TEST_NUM_COMMS *
                                   T8_CMESH_MAX_NUM_OF_TREES))
    % (T8_CMESH_TEST_NUM_COMMS * T8_CMESH_MAX_NUM_OF_TREES);

  t8_debugf
    ("Creating new bigmesh cmesh. eclass=%s,num_trees=%i, comm=%s  \n",
     t8_eclass_to_string[eclass], num_trees, t8_comm_string_list[comm_num]);

  return t8_cmesh_new_bigmesh (eclass, num_trees, comm);
}

/* The function t8_test_create_new_prism_cake_cmesh (int cmesh_id) returns a new_prism_cake cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. The minimal number of trees is 3. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_prism_cake cmesh.
 */
t8_cmesh_t
t8_test_create_new_prism_cake_cmesh (int cmesh_id)
{
 /*
  *    comm    |  num_of_prisms    
  *     0      |       3     
  *     0      |       4     
  *     0      |      ...     
  *     0      |      100     
  *     1      |       3    
  *     1      |       4    
  *    ...     |      ...    
  *     1      |      100  
  */
  const int           min_num_of_prisms = 3;
  int                 comm_num = (cmesh_id / T8_CMESH_MAX_NUM_OF_PRISMS)
    % T8_CMESH_TEST_NUM_COMMS;
  sc_MPI_Comm         comm = t8_comm_list[comm_num];
  int                 num_of_prisms = min_num_of_prisms
    + (cmesh_id % T8_CMESH_MAX_NUM_OF_PRISMS);

  t8_debugf ("Creating new prism cake cmesh. comm=%s, num_of_prisms=%i \n",
             t8_comm_string_list[comm_num], num_of_prisms);

  return t8_cmesh_new_prism_cake (comm, num_of_prisms);
}

/* The function t8_test_create_new_disjoint_bricks_cmesh (int cmesh_id) returns a new_disjoint_bricks cmesh with a unique input for every given id. 
 * The comm is taken from the t8_comm_list. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique new_disjoint_bricks cmesh.
 */
t8_cmesh_t
t8_test_create_new_disjoint_bricks_cmesh (int cmesh_id)
{
  int                 comm_num = cmesh_id % T8_CMESH_TEST_NUM_COMMS;

  sc_MPI_Comm         comm = t8_comm_list[comm_num];

  int                 z_periodic = (cmesh_id / T8_CMESH_TEST_NUM_COMMS)
    % T8_CMESH_BINARY;

  int                 y_periodic = (cmesh_id
                                    / (T8_CMESH_TEST_NUM_COMMS *
                                       T8_CMESH_BINARY))
    % T8_CMESH_BINARY;

  int                 x_periodic = (cmesh_id
                                    / (T8_CMESH_TEST_NUM_COMMS *
                                       T8_CMESH_BINARY * T8_CMESH_BINARY))
    % T8_CMESH_BINARY;

  t8_gloidx_t         num_z = (cmesh_id
                               / (T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY
                                  * T8_CMESH_BINARY * T8_CMESH_BINARY))
    % T8_CMESH_MAX_NUM_XYZ_TREES;

  t8_gloidx_t         num_y = (cmesh_id
                               / (T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY
                                  * T8_CMESH_BINARY * T8_CMESH_BINARY *
                                  T8_CMESH_MAX_NUM_XYZ_TREES))
    % T8_CMESH_MAX_NUM_XYZ_TREES;

  t8_gloidx_t         num_x = (cmesh_id
                               / (T8_CMESH_TEST_NUM_COMMS * T8_CMESH_BINARY
                                  * T8_CMESH_BINARY * T8_CMESH_BINARY
                                  * T8_CMESH_MAX_NUM_XYZ_TREES *
                                  T8_CMESH_MAX_NUM_XYZ_TREES))
    % T8_CMESH_MAX_NUM_XYZ_TREES;

  t8_debugf
    ("Creating new disjoint bricks cmesh. num_x=%li, num_y=%li , num_z=%li , x_periodic=%i, y_periodic=%i, z_periodic=%i, comm=%s \n",
     num_x, num_y, num_z, x_periodic, y_periodic, z_periodic,
     t8_comm_string_list[comm_num]);

  return t8_cmesh_new_disjoint_bricks (num_x, num_y, num_z, x_periodic,
                                       y_periodic, z_periodic, comm);
}

/* The function t8_test_create_cmesh (int cmesh_id) combines all t8_test_create_*_cmesh functions 
 * so that depending on the range the id is in, we get another cmesh type by calling its 
 * t8_test_create_*_cmesh function. 
 * \param [in] cmesh_id The cmesh_id is used to create a unique cmesh.
 */
t8_cmesh_t
t8_test_create_cmesh (int cmesh_id)
{
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_comm_only_cmesh_testcases ()) {
    return t8_test_create_comm_only_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_comm_only_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_hypercube_cmesh_testcases ()) {
    return t8_test_create_new_hypercube_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_hypercube_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_empty_cmesh_testcases ()) {
    return t8_test_create_new_empty_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_empty_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_from_class_cmesh_testcases ()) {
    return t8_test_create_new_from_class_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_from_class_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id <
      t8_get_number_of_new_hypercube_hybrid_cmesh_testcases ()) {
    return t8_test_create_new_hypercube_hybrid_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_hypercube_hybrid_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_periodic_cmesh_testcases ()) {
    return t8_test_create_new_periodic_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_periodic_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_bigmesh_cmesh_testcases ()) {
    return t8_test_create_new_bigmesh_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_bigmesh_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_prism_cake_cmesh_testcases ()) {
    return t8_test_create_new_prism_cake_cmesh (cmesh_id);
  }
  cmesh_id -= t8_get_number_of_new_prism_cake_cmesh_testcases ();
  if (0 <= cmesh_id
      && cmesh_id < t8_get_number_of_new_disjoint_bricks_cmesh_testcases ()) {
    return t8_test_create_new_disjoint_bricks_cmesh (cmesh_id);
  }
  return t8_test_create_comm_only_cmesh (cmesh_id);
}
