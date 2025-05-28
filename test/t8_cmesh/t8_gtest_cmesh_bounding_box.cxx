/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/**
 * \file t8_gtest_cmesh_bounding_box.cxx
 * 
 * Test the computation of the bounding box of a t8_cmesh.
 */

#include <gtest/gtest.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

/**
 * \brief Test fixture for testing the bounding box of a t8_cmesh. Computes a cmesh inside
 * the unit cube, computes the bounding box and checks that it is correct.
 */
class t8_cmesh_bounding_box: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_eclass eclass;
};

TEST_P (t8_cmesh_bounding_box, test_box)
{
  double bounds[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounds);

  const int dim = t8_eclass_to_dimension[eclass];

  for (int idim = 0; idim < dim; ++idim) {
    EXPECT_DOUBLE_EQ (bounds[idim * 2], 0);
    EXPECT_DOUBLE_EQ (bounds[idim * 2 + 1], 1);
  }
}

// Parameterized test cases
INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_bounding_box, t8_cmesh_bounding_box, AllEclasses, print_eclass);
