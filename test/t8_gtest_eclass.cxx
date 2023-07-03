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

/* *INDENT-OFF* */
class t8_gtest_eclass : public testing::TestWithParam<int>{
protected:
  void SetUP() override {
    ieclass = GetParam();
  }
  int ieclass;
};
/* *INDENT-ON* */

TEST (t8_gtest_eclass, eclassCountIs8)
{
  EXPECT_EQ (T8_ECLASS_COUNT, 8);
}

TEST (t8_gtest_eclass, invalid_class)
{
  EXPECT_FALSE (t8_eclass_is_valid ((t8_eclass_t) T8_ECLASS_INVALID));
}

TEST_P (t8_gtest_eclass, dimension)
{
  int                 eclass_dims[8] = { 0, 1, 2, 2, 3, 3, 3, 3 };
  EXPECT_EQ (t8_eclass_to_dimension[ieclass], eclass_dims[ieclass]);
}

TEST_P (t8_gtest_eclass, valid_class)
{
  int                 eclass;
  for (eclass = T8_ECLASS_ZERO; eclass < T8_ECLASS_COUNT; ++eclass) {
    EXPECT_TRUE (t8_eclass_is_valid ((t8_eclass_t) eclass));
  }
}

TEST (t8_gtest_eclass, compare)
{
  int                 eci, ecj;
  for (eci = T8_ECLASS_ZERO; eci < T8_ECLASS_COUNT; ++eci) {
    for (ecj = T8_ECLASS_ZERO; ecj < T8_ECLASS_COUNT; ++ecj) {
      if (eci == ecj) {
        EXPECT_FALSE (t8_eclass_compare
                      ((t8_eclass_t) eci, (t8_eclass_t) ecj));
      }
      else if (t8_eclass_to_dimension[eci] == t8_eclass_to_dimension[ecj]) {
        EXPECT_TRUE (t8_eclass_compare
                     ((t8_eclass_t) eci, (t8_eclass_t) ecj));
      }
    }
  }
}

TEST (t8_gtest_eclass, t8_to_vtk_corner_numbers)
{
  for (int ieclass = T8_ECLASS_ZERO; ieclass < T8_ECLASS_COUNT; ieclass++) {
    const int           num_vertices = t8_eclass_num_vertices[ieclass];
    for (int ivertex = 0; ivertex < num_vertices; ivertex++) {
      const int           vtk_corner_number =
        t8_eclass_t8_to_vtk_corner_number[ieclass][ivertex];
      const int           t8_corner_number =
        t8_eclass_vtk_to_t8_corner_number[ieclass][vtk_corner_number];
      EXPECT_EQ (ivertex, t8_corner_number);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_eclass,
                          testing::Range (T8_ECLASS_ZERO, T8_ECLASS_COUNT));
