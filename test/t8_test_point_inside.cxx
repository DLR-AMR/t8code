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

#include <sc_functions.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_element_cxx.hxx>

/* This function creates a single element of the specified element class.
 * It the creates a bunch of points, some of which lie whithin the element, some
 * not. For each point we call t8_forest_element_point_inside and check whether
 * it gives the correct result. */
/*
 * TODO: - Use more than one refinement level.
 *       - Use barycentric coordinates to create the points, this
 *         spares us to manually check whether a point is inside or not.
 *         (0 <= x_i <= 1 and sum x_i = 1    <=> Point is inside)
 *       - Does the new barycentric coordinate test work with HEX and PRISM?
 */
static void
t8_test_point_inside_level0 (sc_MPI_Comm comm, t8_eclass_t eclass)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_eclass_scheme_c *eclass_scheme;
  t8_element_t       *element;
  double              test_point[3];
  int                 ipoint, icoord;   /* loop variables */
  const int           num_points_per_dim = 5;   /* we construct num_points_per_dim^3 many test points */
  const double        offset = 2.2 / (num_points_per_dim - 1);  /* used to calculate the coordinates of the test points */
  int                 point_is_recognized_as_inside;
  int                 num_corners, icorner;
  double             *tree_vertices;
  double              element_vertices[T8_ECLASS_MAX_CORNERS][3];
  double             *barycentric_coordinates;
  const double        barycentric_range_lower_bound = 0.001;    /* Must be > 0 */
  const double        barycentric_range_upper_bound = 1.1;      /* Should be > 1 */
  int                 num_steps;
  double              step;
  int                 num_points;
  const double        tolerance = 1e-12;        /* Numerical tolerance that we allow in the point inside check */

  default_scheme = t8_scheme_new_default_cxx ();
  /* Construct a cube coarse mesh */
  cmesh = t8_cmesh_new_from_class (eclass, comm);
  /* Build a uniform level 0 forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 0, comm);

  if (t8_forest_get_num_element (forest) > 0) { /* Skip empty forests (occur when executed in parallel) */

    /* Get a pointer to the single element */
    element = t8_forest_get_element (forest, 0, NULL);
    /* Get the vertices of the tree */
    tree_vertices = t8_forest_get_tree_vertices (forest, 0);

    /* Get the associated eclass scheme */
    eclass_scheme = t8_forest_get_eclass_scheme (forest, eclass);

    /* Compute the corner coordinates of the element */
    num_corners = eclass_scheme->t8_element_num_corners (element);
    T8_ASSERT (0 <= num_corners && num_corners <= T8_ECLASS_MAX_CORNERS);       /* Everything else is impossible */
    /* For each corner get its coordinates */
    for (icorner = 0; icorner < num_corners; ++icorner) {
      t8_forest_element_coordinate (forest, 0, element, tree_vertices,
                                    icorner, element_vertices[icorner]);
    }

    /* Allocate the barycentric coordinates */
    barycentric_coordinates = T8_ALLOC (double, num_corners);

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
     * The more corners the smaller the number of stepes in order to limit the
     * number of points that we generate.
     */
    num_points = 10000;         /* Desired number of points */
    /* Compute the number of steps needed */
    num_steps = pow (num_points, 1. / (num_corners - 1));
    /* Corrected number of points due to possible rounding errors in pow */
    num_points = sc_intpow (num_steps, num_corners - 1);
    step =
      (barycentric_range_upper_bound -
       barycentric_range_lower_bound) / (num_steps - 1);
    t8_debugf ("step size %g, steps %i, points %i (corners %i)\n", step,
               num_steps, num_points, num_corners);
    int                 num_in = 0;
    double              dampening;
    for (ipoint = 0; ipoint < num_points; ++ipoint) {
      double              Sum = 0;
      int                 point_is_inside = 1;

      /* Set the coordinates of the test point to 0 */
      for (icoord = 0; icoord < 3; ++icoord) {
        test_point[icoord] = 0;
      }
      for (icorner = 0; icorner < num_corners - 1; ++icorner) {
        int                 this_step =
          (ipoint / sc_intpow (num_steps, icorner)) % num_steps;
        barycentric_coordinates[icorner] =
          barycentric_range_lower_bound + this_step * step;
        dampening = (1 - Sum) * (1 - Sum);
        barycentric_coordinates[icorner] *= dampening;
        Sum += barycentric_coordinates[icorner];

        /* Construct the actual test point */
        for (icoord = 0; icoord < 3; ++icoord) {
          test_point[icoord] +=
            barycentric_coordinates[icorner] *
            element_vertices[icorner][icoord];
        }

        point_is_inside = point_is_inside
          && barycentric_coordinates[icorner] <= 1;
      }
      barycentric_coordinates[num_corners - 1] = 1 - Sum;

      /* Add the last barycentric coordinate to the test point */
      for (icoord = 0; icoord < 3; ++icoord) {
        test_point[icoord] +=
          barycentric_coordinates[num_corners -
                                  1] * element_vertices[num_corners -
                                                        1][icoord];
      }

      point_is_inside = point_is_inside
        && barycentric_coordinates[num_corners - 1] >= 0
        && barycentric_coordinates[num_corners - 1] <= 1;
      num_in += point_is_inside ? 1 : 0;

      /* We now check whether the point inside function correctly sees whether
       * the point is inside the element or not. */
      point_is_recognized_as_inside =
        t8_forest_element_point_inside (forest, 0, element, tree_vertices,
                                        test_point, tolerance);

      SC_CHECK_ABORTF (!point_is_recognized_as_inside == !point_is_inside,
                       "The point (%g,%g,%g) should %s be inside the %s element, but isn't detected as such.",
                       test_point[0], test_point[1], test_point[2],
                       point_is_inside ? "" : "not",
                       t8_eclass_to_string[eclass]);

    }
    t8_debugf ("%i (%.2f%%) of test points are inside the element\n", num_in,
               (100.0 * num_in) / num_points);

    T8_FREE (barycentric_coordinates);
  }                             /* Skip empty forest if */

  t8_forest_unref (&forest);
}

/* In this test we define a triangle in the x-y plane
 * and a point that lies in a triangle that is parallel
 * to this triangle on the z-axis.
 * The point must be correctly identified as lying outside
 * of the triangle.
 */
static void
t8_test_point_inside_specific_triangle ()
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_element_t       *element;
  double              vertices[9] = {
    0., 0., 0.,
    1., 0., 0.,
    1., 1., 0.
  };
  double              test_point[3] = {
    0.3, 0.3, 1
  };
  int                 point_is_inside;
  double             *tree_vertices;
  const double        tolerance = 1e-12;        /* Numerical tolerance that we allow for the point inside check */

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 3);

  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0,
                           sc_MPI_COMM_WORLD);

  if (t8_forest_get_num_element (forest) <= 0) {        /* Skip empty forests (occur when executed in parallel) */
    t8_forest_unref (&forest);
    return;
  }

  element = t8_forest_get_element (forest, 0, NULL);

  /* Get the vertices of the tree */
  tree_vertices = t8_forest_get_tree_vertices (forest, 0);

  point_is_inside =
    t8_forest_element_point_inside (forest, 0,
                                    element, tree_vertices, test_point,
                                    tolerance);

  SC_CHECK_ABORT (!point_is_inside,
                  "The point is wrongly detected as inside the triangle.");

  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;
  int                 ieclass;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf
    ("Testing one specific point with one specific triangle.\n");
  t8_test_point_inside_specific_triangle ();
  for (ieclass = T8_ECLASS_LINE; ieclass < T8_ECLASS_COUNT; ieclass++) {
    if (ieclass != T8_ECLASS_PYRAMID && ieclass != T8_ECLASS_QUAD) {
      /* TODO: - does not work with pyramids yet, since pyramid elements are not implemented.
       *         the point check should work with pyramids, once they are implemented.
       *       - point inside check does not work with quads yet. */
      t8_global_productionf
        ("Testing point finding with eclass %s\n",
         t8_eclass_to_string[ieclass]);
      t8_test_point_inside_level0 (mpic, (t8_eclass_t) ieclass);
    }
  }

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
