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
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_schemes.hxx>

/* This test is executed on a subcommunicator of exactly 2 procs, because it demonstrates a configuration that is currently not working. See https://github.com/DLR-AMR/t8code/issues/825.
 * A partitioned square of uniform refinement level 1 is adapted once, where only the lower half is refined.
 * Then the mesh is adapted again, where only the upper half of the lower elements is deleted.
 * Now, no element actually has a face connection to an existing element on the other proc, thus the Ghost structure should be empty.
 * Within the bug discussed in https://github.com/DLR-AMR/t8code/issues/825 they send a message anyway.
 * Once, https://github.com/DLR-AMR/t8code/issues/825 is resolved, no messages should be send and the test should pass.
 */

/* refine elements, whose lower left y coordinate is 0, so that the lowest row is refined
 * delete elements, whose lower left y coordniate is 0.25, so that the second row is deleted
 */
static int
test_adapt_holes (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, const t8_eclass_t tree_class,
                  t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family, const int num_elements,
                  t8_element_t *elements[])
{
  double coordinates[3];
  t8_forest_element_coordinate (forest_from, which_tree, elements[0], 0, coordinates);
  if (fabs (coordinates[1]) < T8_PRECISION_EPS)
    return 1;  // refine the lowest row
  if (fabs (coordinates[1] - 0.25) < T8_PRECISION_EPS)
    return -2;  // delete the higher row in the lower quadrants
  return 0;
}

class DISABLED_forest_ghost_exchange_holes: public testing::TestWithParam<std::tuple<int, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    /* adjust communicator size, we split the comm here, because we know that the test fails on 2 procs*/
    int size, rank, color, key;

    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
    sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);

    if (rank < 2) {
      color = 0;
      key = rank;
    }
    else {
      color = sc_MPI_UNDEFINED;
      key = -1;
    }
    sc_MPI_Comm_split (sc_MPI_COMM_WORLD, color, key, &comm);
    T8_ASSERT (rank < 2 || comm == sc_MPI_COMM_NULL);

    if (comm != sc_MPI_COMM_NULL) {
      sc_MPI_Comm_size (comm, &size);
      T8_ASSERT (size <= 2);
      const int scheme_id = std::get<0> (GetParam ());
      scheme = create_from_scheme_id (scheme_id);
      eclass = std::get<1> (GetParam ());
      /* Construct a cmesh */
      cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    }
    else {
      T8_ASSERT (rank >= 2);
    }
  }
  void
  TearDown () override
  {
    //cmesh and scheme are freed taken and freed by forest
    if (comm != sc_MPI_COMM_NULL) {
      sc_MPI_Comm_free (&comm);
    }
    sc_MPI_Barrier (sc_MPI_COMM_WORLD);
    t8_cmesh_unref (&cmesh);
    scheme->unref ();
  }
  sc_MPI_Comm comm;
  const t8_scheme *scheme;
  t8_cmesh_t cmesh;
  t8_eclass_t eclass;
};

TEST_P (DISABLED_forest_ghost_exchange_holes, errorTest)
{
  /* This test tests the functionality described in Issue: https://github.com/DLR-AMR/t8code/issues/825
  */
  if (comm != sc_MPI_COMM_NULL) {
    const int level = 1;
    const int execute_ghost = 1;
    t8_cmesh_ref (cmesh);
    scheme->ref ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 1, comm);
    forest = t8_forest_new_adapt (forest, test_adapt_holes, 0, execute_ghost, NULL);
    forest = t8_forest_new_adapt (forest, test_adapt_holes, 0, execute_ghost, NULL);
    t8_forest_unref (&forest);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_delete, DISABLED_forest_ghost_exchange_holes, AllSchemes);
