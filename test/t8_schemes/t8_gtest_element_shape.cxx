/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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
#include <t8_eclass/t8_eclass.h>
#include <t8_element/t8_element_shape.h>
#include <test/t8_gtest_macros.hxx>

struct gtest_shapes: public testing::TestWithParam<t8_eclass_t>
{
 protected:
  void
  SetUp () override
  {
    shape = (t8_element_shape_t) GetParam ();
  }
  t8_element_shape_t shape;
};

TEST_P (gtest_shapes, t8_to_vtk_corner_number)
{
  const int corner_nums = t8_element_shape_num_vertices (shape);
  for (int icorner = 0; icorner < corner_nums; icorner++) {
    const int vtk_corner_number = t8_element_shape_t8_to_vtk_corner_number (shape, icorner);
    const int control_corner = t8_element_shape_vtk_to_t8_corner_number (shape, vtk_corner_number);
    EXPECT_EQ (icorner, control_corner);
  }
}

TEST_P (gtest_shapes, shape_num_faces)
{
  const int dim = t8_eclass_to_dimension[(t8_eclass_t) shape];
  EXPECT_LE (t8_element_shape_num_faces (shape), t8_element_shape_max_num_faces (dim));
}

TEST_P (gtest_shapes, shape_compare)
{
  for (int ishape = T8_ECLASS_ZERO; ishape < T8_ECLASS_COUNT; ++ishape) {
    if (ishape == shape) {
      EXPECT_FALSE (t8_element_shape_compare ((t8_element_shape_t) ishape, shape));
    }
    else if (t8_eclass_to_dimension[ishape] == t8_eclass_to_dimension[shape]) {
      EXPECT_TRUE (t8_element_shape_compare ((t8_element_shape_t) ishape, shape));
    }
  }
}

TEST_P (gtest_shapes, shape_string)
{
  const char* shape_string = t8_element_shape_to_string (shape);
  EXPECT_TRUE (shape_string != nullptr);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_shapes, gtest_shapes, AllEclasses, print_eclass);
