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
#include <sc_functions.h>
#include <t8_eclass.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>

/* In this test we define a triangle in the x-y plane
 * and a point that lies in a triangle that is parallel
 * to this triangle on the z-axis.
 * The point must be correctly identified as lying outside
 * of the triangle.
 */
TEST (t8_point_inside, test_point_inside_specific_triangle)
{
  t8_cmesh_t cmesh;

  /* clang-format off */
  double vertices[9] = { 0., 0., 0., 
                         1., 0., 0., 
                         1., 1., 0. };
  /* clang-format on */
  double test_point[3] = { 0.3, 0.3, 1 };
  const double tolerance = 1e-12; /* Numerical tolerance that we allow for the point inside check */

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  /* We use standard linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, sc_MPI_COMM_WORLD);

  if (t8_forest_get_local_num_elements (forest) <= 0) {
    /* Skip empty forests (can occur when executed in parallel) */
    t8_forest_unref (&forest);
    GTEST_SKIP ();
  }

  t8_element_t *element = t8_forest_get_element (forest, 0, NULL);

  int point_is_inside;
  t8_forest_element_points_inside (forest, 0, element, test_point, 1, &point_is_inside, tolerance);
  ASSERT_FALSE (point_is_inside) << "The point is wrongly detected as inside the triangle.";
  t8_forest_unref (&forest);
}

/* In this test we define a quad in the x-y plane
 * and a point that lies in a quad that is parallel
 * to the first quad on the z-axis.
 * The point must be correctly identified as lying outside
 * of the quad.
 */
TEST (t8_point_inside, test_point_inside_specific_quad)
{
  t8_cmesh_t cmesh;

  /* clang-format off */
  double vertices[12] = { 0., 0., 0., 
                          1., 0., 0., 
                          0., 1., 0., 
                          1., 1., 0. };
  /* clang-format on */
  double test_point[3] = { 0.3, 0.3, 1 };
  const double tolerance = 1e-12; /* Numerical tolerance that we allow for the point inside check */

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  /* We use standard linear geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 0, 0, sc_MPI_COMM_WORLD);

  if (t8_forest_get_local_num_elements (forest) <= 0) {
    /* Skip empty forests (can occur when executed in parallel) */
    t8_forest_unref (&forest);
    GTEST_SKIP ();
  }

  t8_element_t *element = t8_forest_get_element (forest, 0, NULL);

  int point_is_inside;
  t8_forest_element_points_inside (forest, 0, element, test_point, 1, &point_is_inside, tolerance);

  ASSERT_FALSE (point_is_inside) << "The point is wrongly detected as inside the quad.";

  t8_forest_unref (&forest);
}

class geometry_point_inside: public testing::TestWithParam<std::tuple<t8_eclass, int, int>> {
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    level = std::get<1> (GetParam ());
    use_axis_aligned_geom = std::get<2> (GetParam ());

    /* Construct a cube coarse mesh */
    if (use_axis_aligned_geom && (eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX)) {
      /* clang-format off */
      const double boundaries[24] = { 
        0, 0, 0, 
        1, 0, 0, 
        0, 1, 0, 
        1, 1, 0, 
        0, 0, 1, 
        1, 0, 1, 
        0, 1, 1,
        1, 1, 1 
      };
      /* clang-format on */
      cmesh = t8_cmesh_new_hypercube_pad (eclass, sc_MPI_COMM_WORLD, boundaries, 1, 1, 1, use_axis_aligned_geom);
    }
    else {
      cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
    }
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
  int level;
  int use_axis_aligned_geom;
  t8_cmesh_t cmesh;
};

TEST_P (geometry_point_inside, test_point_inside)
{
  double element_vertices[T8_ECLASS_MAX_CORNERS][3];
  const double barycentric_range_lower_bound = 0.1; /* Must be > 0 */
  const double barycentric_range_upper_bound = 1.1; /* Should be > 1 */
  /* Numerical tolerance that we allow in the point inside check */
  const double tolerance = 1e-12;
  const int num_points_to_generate = 64; /* Desired number of test points per element. 
                                                         * The actual number is computed to be close to this and
                                                         * will always be >= 2^num_corners (thus >= 128 for hex). */

  t8_debugf ("Testing eclass %s, uniform level %i with approx. %i points per element.\n", t8_eclass_to_string[eclass],
             level, num_points_to_generate);
  const t8_scheme *default_scheme = t8_scheme_new_default ();

  /* We translate the coordinates of the cmesh to create a non-standard case.
   * In particular, we want the 1D and 2D elements to move outside of axis
   * perpendicular planes. We do this to detect possible errors in the point
   * finding function that may be biased towards certrain coordinate axis.
   */
  double *const tree_vertices = t8_cmesh_get_tree_vertices (cmesh, 0);
  const int num_vertices = (use_axis_aligned_geom && (eclass == T8_ECLASS_HEX || eclass == T8_ECLASS_QUAD))
                             ? 2
                             : t8_eclass_num_vertices[eclass];
  /* Translate all points by the same vector to move the element a bit. */
  double translate_all_points[3] = { -0.1, 0.3, 0.15 };
  t8_cmesh_translate_coordinates (tree_vertices, tree_vertices, num_vertices, translate_all_points);
  /* Translate points 0 and 1 (if it exists) extra in order to move the 2D elements
   * and 3D faces outside of axis perpendicular planes. */
  if (!use_axis_aligned_geom) {
    double translate_points_0_1[3] = { 0.1, -0.1, 0.3 };
    t8_cmesh_translate_coordinates (tree_vertices, tree_vertices, 1, translate_points_0_1);
    if (num_vertices > 2) {
      t8_cmesh_translate_coordinates (tree_vertices + 3, tree_vertices + 3, 1, translate_points_0_1);
    }
  }

  /* Build a uniform forest */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, sc_MPI_COMM_WORLD);

  const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    t8_log_indent_push ();
    const t8_locidx_t num_elements = t8_forest_get_tree_num_elements (forest, itree);
    /* Get the associated eclass scheme */
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    for (t8_locidx_t ielement = 0; ielement < num_elements; ++ielement) {
      /* Get a pointer to the element */
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);

      /* Compute the corner coordinates of the element */
      const int num_corners = scheme->element_get_num_corners (tree_class, element);
      /* For each corner get its coordinates */
      for (int icorner = 0; icorner < num_corners; ++icorner) {
        t8_forest_element_coordinate (forest, itree, element, icorner, element_vertices[icorner]);
      }

      /* Allocate the barycentric coordinates */
      double *barycentric_coordinates = T8_ALLOC (double, num_corners);

      /* Fill the barycentric coordinates with test values */
      /* 1-Sum 0 0 0  ... start
       * 1-Sum 0 0 0  ... start + step
       *   .
       *   .
       * 1-Sum 0 0 0  ... end
       * 1-Sum 0 0 0    start start
       * 1-Sum 0 0 0    start + step start
       *  .
       *  .
       *  .
       *
       * Thus, we have num_steps^(num_corners - 1) many points.
       * if 0 <= i < num_points then the j-th coordinate has to be
       * set to start + (i / num_steps^j) % num_steps * step
       *
       * We multiply this entry with a dampening factor that is smaller the smaller
       * the sum of the already set entries is.
       * This dampening increases the percentage of points that lie inside the element.
       *
       * The number of steps depends on the number of corners.
       * The more corners the smaller the number of steps in order to limit the
       * number of points that we generate.
       */
      /* Compute the number of steps needed */
      int num_steps = pow (num_points_to_generate, 1. / (num_corners - 1));
      if (num_steps <= 1) {
        /* Correct to at least 2 steps in order to prevent division by 0 */
        num_steps = 2;
      }
      /* Corrected number of points due to possible rounding errors in pow */
      const int min_points_outside = 6;
      const int num_points = sc_intpow (num_steps, num_corners - 1);
      const int total_points = num_points + min_points_outside;
      double *test_point = T8_ALLOC_ZERO (double, total_points * 3);
      int *point_is_inside = T8_ALLOC (int, total_points);
      int *point_is_recognized_as_inside = T8_ALLOC (int, total_points);
      double step = (barycentric_range_upper_bound - barycentric_range_lower_bound) / (num_steps - 1);
      //t8_debugf ("step size %g, steps %i, points %i (corners %i)\n", step,
      //           num_steps, num_points, num_corners);
      int num_in = 0;
      double dampening;

      /* For each test point, build coordinates and test whether it
       * is correctly detected to be in/outside the element. */
      for (int ipoint = 0; ipoint < num_points; ++ipoint) {
        double Sum = 0;
        point_is_inside[ipoint] = 1;

        /* Set the coordinates of the test point to 0 */
        for (int icoord = 0; icoord < 3; ++icoord) {
          test_point[ipoint * 3 + icoord] = 0;
        }
        for (int icorner = 0; icorner < num_corners - 1; ++icorner) {
          int this_step = (ipoint / sc_intpow (num_steps, icorner)) % num_steps;
          /* Set barycentric coordinates */
          barycentric_coordinates[icorner] = barycentric_range_lower_bound + this_step * step;
          dampening = (1 - Sum) * (1 - Sum);
          barycentric_coordinates[icorner] *= dampening;
          Sum += barycentric_coordinates[icorner];

          /* Construct the actual test point from the barycentric coordinate */
          for (int icoord = 0; icoord < 3; ++icoord) {
            test_point[ipoint * 3 + icoord] += barycentric_coordinates[icorner] * element_vertices[icorner][icoord];
          }

          /* The point is inside if and only if all barycentric coordinates are >= 0. */
          point_is_inside[ipoint] = point_is_inside[ipoint] && barycentric_coordinates[icorner] >= 0;
        }
        /* Ensure that sum over all bar. coordinates is 1 */
        barycentric_coordinates[num_corners - 1] = 1 - Sum;

        /* Add this last barycentric coordinate to the test point */
        for (int icoord = 0; icoord < 3; ++icoord) {
          test_point[ipoint * 3 + icoord]
            += barycentric_coordinates[num_corners - 1] * element_vertices[num_corners - 1][icoord];
        }

        /* The point is inside if and only if all barycentric coordinates are >= 0. */
        point_is_inside[ipoint] = point_is_inside[ipoint] && barycentric_coordinates[num_corners - 1] >= 0;

        num_in += point_is_inside ? 1 : 0;
      }
      /* Create a set of points that are outside of the cube. 
       * The coordinates are given by the corner points except for one coordinate. For each side of the cube
       * we place one point outside of it. */
      for (int side = 0; side < min_points_outside; ++side) {
        test_point[3 * (num_points + side)] = (side == 0) ? -2.0 : (side == 1) ? 2.0 : -0.1;
        test_point[3 * (num_points + side) + 1] = (side == 2) ? -2.0 : (side == 3) ? 2.0 : 0.3;
        test_point[3 * (num_points + side) + 2] = (side == 4) ? -2.0 : (side == 5) ? 2.0 : 0.15;
        point_is_inside[num_points + side] = false;
      }
      /* We now check whether the point inside function correctly sees whether
         * the point is inside the element or not. */
      t8_forest_element_points_inside (forest, 0, element, test_point, total_points, point_is_recognized_as_inside,
                                       tolerance);
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        ASSERT_EQ (point_is_recognized_as_inside[ipoint], point_is_inside[ipoint])
          << "Testing point #" << ipoint << "(" << test_point[0] << "," << test_point[1] << "," << test_point[2]
          << ") should " << (point_is_inside[ipoint] ? "" : "not ") << "be inside the " << t8_eclass_to_string[eclass]
          << " element, but is not detected as such.";
      } /* End loop over points. */
      t8_debugf ("%i (%.2f%%) test points are inside the element\n", num_in, (100.0 * num_in) / num_points);
      T8_FREE (barycentric_coordinates);
      T8_FREE (test_point);
      T8_FREE (point_is_inside);
      T8_FREE (point_is_recognized_as_inside);
    } /* End loop over elements */
  }   /* End loop over trees */
  t8_forest_unref (&forest);
  t8_log_indent_pop ();
}

auto print_test = [] (const testing::TestParamInfo<std::tuple<t8_eclass, int, int>> &info) {
  const t8_eclass_t eclass = std::get<0> (info.param);
  const int level = std::get<1> (info.param);
  const int use_axis_aligned_geom = std::get<2> (info.param);

  const std::string geom = use_axis_aligned_geom ? std::string ("AxisAligned") : std::string ("Linear");

  std::string name
    = std::string (t8_eclass_to_string[eclass]) + std::string ("_") + std::to_string (level) + std::string ("_") + geom;
  return name;
};

#if T8CODE_TEST_LEVEL >= 1
INSTANTIATE_TEST_SUITE_P (t8_gtest_point_inside, geometry_point_inside,
                          testing::Combine (testing::Range (T8_ECLASS_LINE, T8_ECLASS_QUAD), testing::Range (0, 4),
                                            testing::Range (0, 2)),
                          print_test);

#else
INSTANTIATE_TEST_SUITE_P (t8_gtest_point_inside, geometry_point_inside,
                          testing::Combine (testing::Range (T8_ECLASS_LINE, T8_ECLASS_COUNT), testing::Range (0, 6),
                                            testing::Range (0, 2)),
                          print_test);
#endif
