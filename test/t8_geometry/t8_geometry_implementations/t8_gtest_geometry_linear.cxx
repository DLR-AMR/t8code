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

/** \file t8_gtest_geometry_linear.cxx
 * This file contains tests for the linear geometry of t8code.
 */

#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_element.h>
#include <t8_types/t8_vec.hxx>

class geometry_test: public testing::TestWithParam<std::tuple<int, t8_eclass>> {
 public:
  static void
  SetUpTestSuite ()
  {
    seed = time (NULL); /* RNG seed */
  }
  static int seed;

 protected:
  void
  SetUp () override
  {
    const int geom_int = std::get<0> (GetParam ());
    eclass = std::get<1> (GetParam ());
    t8_cmesh_init (&cmesh);
    if (geom_int == T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED
        && !(eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX)) {
      GTEST_SKIP ();
    }

    const int num_vertices = t8_eclass_num_vertices[eclass];
    t8_cmesh_set_tree_class (cmesh, 0, eclass);
    double *vertices = T8_ALLOC_ZERO (double, num_vertices *T8_ECLASS_MAX_DIM);
    for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
      for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
        vertices[i_vertex * T8_ECLASS_MAX_DIM + dim] = t8_element_corner_ref_coords[eclass][i_vertex][dim];
      }
    }
    switch (geom_int) {
    case T8_GEOMETRY_TYPE_LINEAR:
      geom = new t8_geometry_linear ();
      t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
      break;
    case T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED:
      geom = new t8_geometry_linear_axis_aligned ();
      /* Copy last vertex to the second position*/
      vertices[3] = vertices[3 * (num_vertices - 1)];
      vertices[4] = vertices[3 * (num_vertices - 1) + 1];
      vertices[5] = vertices[3 * (num_vertices - 1) + 2];

      t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (cmesh);
      break;
    default:
      break;
    }
    t8_cmesh_set_tree_vertices (cmesh, 0, vertices,
                                geom_int == T8_GEOMETRY_TYPE_LINEAR ? t8_eclass_num_vertices[eclass] : 2);
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
    T8_FREE (vertices);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }
  t8_eclass_t eclass;
  t8_geometry_with_vertices *geom;
  t8_cmesh_t cmesh;
};

int geometry_test::seed;

/* Check whether the linear axis aligned geometry map is correct.
 * We create a cmesh of one tree and unit line, quad or hex
 * geometry. We then create random points in a reference tree and
 * check whether the evaluation is correct. */
TEST_P (geometry_test, cmesh_geometry)
{
  /* TODO: Add a test for the jacobian, as soon as its implemented. */

  /* Create random points in [0,1]^d and check if they are mapped correctly. */

  t8_3D_point point_mapped;
  const t8_geometry_c *cmesh_geom;

  /* Double check that the geometry is the linear axis aligned geometry. */
  cmesh_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_TRUE (cmesh_geom != NULL) << "Could not get cmesh's geometry.";
  ASSERT_EQ (cmesh_geom->t8_geom_get_hash (), geom->t8_geom_get_hash ())
    << "cmesh's geometry is not the expected geometry.";

  srand (seed);
  for (int ipoint = 0; ipoint < T8_NUM_SAMPLE_POINTS; ++ipoint) {
    t8_3D_point point ({ 0.0, 0.0, 0.0 });
    /* Compute random coordinates in [0,1].
     * These are seen as reference coordinates in the single
     * cmesh tree. Our geometry will map them into the physical
     * space. Since this space is also [0,1] and the cmesh only
     * has one tree, the mapped coordinates must be the same as the 
     * reference coordinates. */
    for (int idim = 0; idim < t8_eclass_to_dimension[eclass]; ++idim) {
      point[idim] = (double) rand () / RAND_MAX;
    }

    /* Evaluate the geometry */
    t8_geometry_evaluate (cmesh, 0, point.data (), 1, point_mapped.data ());
    /* Check that the first dim coordinates are the same */
    EXPECT_VEC_EQ (point, point_mapped, T8_PRECISION_SQRT_EPS);
  }
}

auto print_test = [] (const testing::TestParamInfo<std::tuple<int, t8_eclass>> &info) {
  std::string name;
  const int geom_int = std::get<0> (info.param);
  const t8_eclass_t eclass = std::get<1> (info.param);
  switch (geom_int) {
  case T8_GEOMETRY_TYPE_LINEAR:
    name = std::string ("linear_geometry_");
    break;
  case T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED:
    name = std::string ("linear_axis_aligned_geometry_");
    break;
  default:
    name = std::string ("geometry_not_allowed_");
    break;
  }
  name += std::string (t8_eclass_to_string[eclass]);
  return name;
};

INSTANTIATE_TEST_SUITE_P (
  t8_gtest_geometry, geometry_test,
  ::testing::Combine (::testing::Values (T8_GEOMETRY_TYPE_LINEAR, T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED), AllEclasses),
  print_test);

#ifdef T8_ENABLE_DEBUG
TEST (test_geometry_linear, incompatible_geometry)
{
  t8_cmesh_t cmesh;

  t8_debugf ("Testing geometry compatibility checking for linear axis aligned geometry.\n");

  /* Build a simple set geometries for the tree. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_QUAD], 4);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  /* Register the t8_geometry_linear_axis_aligned geometry to this cmesh. */
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (cmesh);
  /* Should return true since the t8_geometry_linear_axis_aligned geometry is compatible with quads. */
  ASSERT_TRUE (t8_cmesh_validate_geometry (cmesh));
  t8_cmesh_destroy (&cmesh);

  /* Build a simple set geometries for the tree. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_TRIANGLE], 3);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 1, *t8_element_corner_ref_coords[T8_ECLASS_QUAD], 4);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  /* Register the linear axis aligned geometry to this cmesh.
   * We register it after committing because it would throw an assertion and do not have death tests.*/
  t8_cmesh_register_geometry<t8_geometry_linear_axis_aligned> (cmesh);
  /* Check validity after committing to circumvent the assertion.
   * Should return false since the t8_geometry_linear_axis_aligned geometry is not compatible with triangles. */
  ASSERT_FALSE (t8_cmesh_validate_geometry (cmesh));
  t8_cmesh_destroy (&cmesh);
}
#endif /* T8_ENABLE_DEBUG */
