/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/*
  This test checks that creating a hypercube coarse mesh on processor 0
  and broadcasting it results in the same coarse mesh as creating it
  on all processes directly
*/

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

class cmesh_hypercube: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();

    cmesh_bcast = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 1, 0, 0);
    cmesh_check = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh_bcast);
    t8_cmesh_unref (&cmesh_check);
  }
  t8_cmesh_t cmesh_bcast, cmesh_check;
  t8_eclass eclass;
};

TEST_P (cmesh_hypercube, bcast_equal_no_bcast)
{
  EXPECT_TRUE (t8_cmesh_is_equal (cmesh_bcast, cmesh_check));
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_bcast, cmesh_hypercube, AllEclasses, print_eclass);
