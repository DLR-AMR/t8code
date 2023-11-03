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
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* This test
 */

static int test_adapt_holes (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                       t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
{
    double coordinates[3];
    double eps = 0.0001;
    t8_forest_element_coordinate (forest_from, which_tree, elements[0], 0, coordinates);
    if (fabs(coordinates[1]) < eps) return 1; // refine the lowest row 
    if (fabs(coordinates[1] - 0.25) < eps) return -2; // refine the lowest row 
    return 0;
}

class forest_ghost_exchange_holes: public testing::Test {
 protected:
  void
  SetUp () override
  {
    /* adjust communicator size */
    int size, rank, color, key;

    sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &size);
    sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
    t8_debugf("size is %i\n", size);

    if (rank < 2){
        color = 0;
        key = rank;
    }else{
        color = sc_MPI_UNDEFINED;
        key = -1;
    }

    t8_debugf("Color: %i, Key: %i\n", color, key);

    sc_MPI_Comm_split(sc_MPI_COMM_WORLD, color, key, &comm);

    if(comm != sc_MPI_COMM_NULL){
        sc_MPI_Comm_size(comm, &size);
        t8_debugf("\nNew size is: %i\n\n", size);
        scheme = t8_scheme_new_default_cxx ();
        /* Construct a cmesh */
        cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    }else{
        t8_debugf("\nThis proc has comm NULL \n\n");
    }

  }
  void
  TearDown () override
  {
    //cmesh and scheme are freed taken and freed by forest
    if(comm != sc_MPI_COMM_NULL){
      sc_MPI_Comm_free(&comm);
    }
    sc_MPI_Barrier(sc_MPI_COMM_WORLD);
  }
  sc_MPI_Comm comm;
  t8_scheme_cxx_t *scheme;
  t8_cmesh_t cmesh;
};

TEST_F (forest_ghost_exchange_holes, errorTest)
{
  if(comm != sc_MPI_COMM_NULL){
    int level = 1;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 1, comm);
    forest = t8_forest_new_adapt (forest, test_adapt_holes, 0, 1, NULL);
    forest = t8_forest_new_adapt (forest, test_adapt_holes, 0, 1, NULL);
    t8_forest_unref(&forest);
  }
}
