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

#include <gtest/gtest.h>
#include <t8_eclass.h>

class gtest_eclass: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    ieclass = GetParam ();
  }
  int ieclass;
};

TEST (gtest_eclass, eclassCountIs8)
{
  EXPECT_EQ (T8_ECLASS_COUNT, 8);
}

TEST (gtest_eclass, invalid_class)
{
  EXPECT_FALSE (t8_eclass_is_valid ((t8_eclass_t) T8_ECLASS_INVALID));
}

TEST_P (gtest_eclass, dimension)
{
  const int eclass_dims[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  EXPECT_EQ (t8_eclass_to_dimension[ieclass], eclass_dims[ieclass]);
}

TEST_P (gtest_eclass, valid_class)
{
  EXPECT_TRUE (t8_eclass_is_valid ((t8_eclass_t) ieclass));
}

TEST_P (gtest_eclass, compare)
{
  for (int jeclass = T8_ECLASS_ZERO; jeclass < T8_ECLASS_COUNT; ++jeclass) {
    if (ieclass == jeclass) {
      EXPECT_FALSE (t8_eclass_compare ((t8_eclass_t) ieclass, (t8_eclass_t) jeclass));
    }
    else if (t8_eclass_to_dimension[ieclass] == t8_eclass_to_dimension[jeclass]) {
      EXPECT_TRUE (t8_eclass_compare ((t8_eclass_t) ieclass, (t8_eclass_t) jeclass));
    }
  }
}

TEST_P (gtest_eclass, t8_to_vtk_corner_numbers)
{
  const int num_vertices = t8_eclass_num_vertices[ieclass];
  for (int ivertex = 0; ivertex < num_vertices; ivertex++) {
    const int vtk_corner_number = t8_eclass_t8_to_vtk_corner_number[ieclass][ivertex];
    const int t8_corner_number = t8_eclass_vtk_to_t8_corner_number[ieclass][vtk_corner_number];
    EXPECT_EQ (ivertex, t8_corner_number);
  }
}

TEST (gtest_eclass, eclass_order)
{
  /* 2D order of eclasses */
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_TRIANGLE, T8_ECLASS_QUAD), -1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE), 1);

  /* 3d order of eclasses */
  /* TET < HEX < PRISM < PYRAMID */
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_TET, T8_ECLASS_HEX), -1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_TET, T8_ECLASS_PRISM), -1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_TET, T8_ECLASS_PYRAMID), -1);

  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_HEX, T8_ECLASS_PRISM), -1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_HEX, T8_ECLASS_PYRAMID), -1);

  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PRISM, T8_ECLASS_PYRAMID), -1);

  /* PYRAMID > PRISM > HEX > TET */
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PYRAMID, T8_ECLASS_PRISM), 1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PYRAMID, T8_ECLASS_HEX), 1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PYRAMID, T8_ECLASS_TET), 1);

  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PRISM, T8_ECLASS_HEX), 1);
  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_PRISM, T8_ECLASS_TET), 1);

  EXPECT_EQ (t8_eclass_compare (T8_ECLASS_HEX, T8_ECLASS_TET), 1);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_eclass, gtest_eclass, testing::Range ((int) T8_ECLASS_ZERO, (int) T8_ECLASS_COUNT));
