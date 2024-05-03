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

/* This gtest checks whether a cmesh of a hypercube is committed and whether for each tree is face consistent in respect to
 * the creation of the cmesh with or without broadcast/ partition. */

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

/* Create class for parameterized Test with multiple test parameters */
class cmesh_hypercube_trees: public testing::TestWithParam<std::tuple<t8_eclass, int, int>> {
 protected:
  /* SetUp the test parameters (eclass, bcast and partition) and define the test value cmesh. */
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    bcast = std::get<1> (GetParam ());
    partition = std::get<2> (GetParam ());

    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, bcast, partition, 0);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_eclass eclass;
  int bcast;
  int partition;
};

/* Test wheater the created cmesh of a hypercube is committed and its are face consistent. Test will fail, if one of these is false. */
TEST_P (cmesh_hypercube_trees, check_cmesh_and_its_trees)
{

  EXPECT_TRUE (t8_cmesh_is_committed (cmesh));
  EXPECT_TRUE (t8_cmesh_trees_is_face_consistent (cmesh, cmesh->trees));
  ASSERT_EQ (t8_cmesh_get_dimension (cmesh), t8_eclass_to_dimension[eclass]) << "Wrong dimension set for cmesh.";
}

/* Use the testing range for eclass with [T8_ECLASS_ZERO, T8_ECLASS_COUNT]. For the generation of the cmesh with or withaout broadcast
 * the booleans 0 and 1 are used. Analogue with partition. */
INSTANTIATE_TEST_SUITE_P (t8_gtest_hypercube, cmesh_hypercube_trees,
                          testing::Combine (AllEclasses, testing::Values (0, 1), testing::Values (0, 1)));
