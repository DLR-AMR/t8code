/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_vec.h>
#include <test/t8_gtest_macros.hxx>

class gtest_eclass: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    ieclass = GetParam ();
  }
  t8_eclass_t ieclass;
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

TEST_P (gtest_eclass, eclass_face_orientation)
{
  const int num_faces = t8_eclass_num_faces[ieclass];
  /* For dim <= 2 all faces are 0 oriented. */
  if (t8_eclass_to_dimension[ieclass] <= 2) {
    for (int iface = 0; iface < num_faces; iface++) {
      EXPECT_EQ (0, t8_eclass_face_orientation[ieclass][iface]);
    }
  }
  else {
    /* Compute the normal of each face and the signed distance between the normal and the vector
     * from vertex 0 to the centroid. If the sign is positive the face is 0 oriented, 
     * ow it is 1 oriented. */
    t8_cmesh_t cmesh = t8_cmesh_new_from_class ((t8_eclass_t) ieclass, sc_MPI_COMM_WORLD);
    const double *vertices = t8_cmesh_get_tree_vertices (cmesh, 0);
    for (int iface = 0; iface < num_faces; iface++) {
      /* clang-format off */
      const int v[3] = {  t8_face_vertex_to_tree_vertex[ieclass][iface][0], 
                          t8_face_vertex_to_tree_vertex[ieclass][iface][1],
                          t8_face_vertex_to_tree_vertex[ieclass][iface][2] };
      /* clang-format on */
      double vec_0[3];
      double vec_1[3];
      t8_vec_axpyz (vertices + 3 * v[1], vertices + 3 * v[0], vec_0, -1.0);
      t8_vec_axpyz (vertices + 3 * v[2], vertices + 3 * v[0], vec_1, -1.0);
      double normal[3];
      t8_vec_cross (vec_0, vec_1, normal);
      double centroid[3] = { 0.0 };
      for (int ivertex = 0; ivertex < t8_eclass_num_vertices[ieclass]; ivertex++) {
        centroid[0] += vertices[3 * ivertex];
        centroid[1] += vertices[3 * ivertex + 1];
        centroid[2] += vertices[3 * ivertex + 2];
      }
      t8_vec_ax (centroid, 1.0 / t8_eclass_num_vertices[ieclass]);
      t8_vec_axpy (vertices + 3 * v[0], centroid, -1.0);
      const double c_n = t8_vec_dot (centroid, normal);
      const int orientation = c_n > 0 ? 0 : 1;
      EXPECT_EQ (orientation, t8_eclass_face_orientation[ieclass][iface]);
    }
    t8_cmesh_destroy (&cmesh);
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

INSTANTIATE_TEST_SUITE_P (t8_gtest_eclass, gtest_eclass, AllEclasses, print_eclass);
