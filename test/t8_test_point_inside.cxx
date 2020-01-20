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

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>

static void
t8_test_point_inside_level0 (sc_MPI_Comm comm, t8_eclass_t eclass)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_element_t       *element;
  double              test_point[3];
  int                 ipoint, jpoint, kpoint;   /* loop variables */
  const int           num_points_per_dim = 5;   /* we construct num_points_per_dim^3 many test points */
  const double        offset = 2.2 / (num_points_per_dim - 1);  /* used to calculate the coordinates of the test points */
  int                 point_is_inside;
  double             *tree_vertices;

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

    /* Create a bunch of points and check whether they are inside the element */
    for (ipoint = 0; ipoint < num_points_per_dim; ++ipoint) {
      test_point[0] = 1.101 - ipoint * offset;  /* we add 0.001 in order to avoid boundary cases where the
                                                   inside check may give a different result than our own check
                                                   du to rounding errors. */
      for (jpoint = 0; jpoint < num_points_per_dim; ++jpoint) {
        if (t8_eclass_to_dimension[eclass] >= 2) {
          /* Only set the y coordinate to non-zero if we test elements of dimension 2 or 3 */
          test_point[1] = 1.102 - jpoint * offset;
        }
        else {
          test_point[1] = 0.0;
        }
        for (kpoint = 0; kpoint < num_points_per_dim; ++kpoint) {
          if (t8_eclass_to_dimension[eclass] == 3) {
            /* only set the z coordinate to non-zero if we test elements of dimension 3 */
            test_point[2] = 1.103 - kpoint * offset;
          }
          else {
            test_point[2] = 0.0;
          }
          point_is_inside =
            t8_forest_element_point_inside (forest, 0, element, tree_vertices,
                                            test_point);

          /* We now manually check whether the point is inside the element or not. */
          switch (eclass) {
          case T8_ECLASS_QUAD:
            /* The point is inside if all its x and y coordinates are in [0, 1] and z = 0 */
            if (test_point[0] >= 0 && test_point[1] >= 0 && test_point[2] == 0
                && test_point[0] <= 1 && test_point[1] <= 1
                && test_point[2] == 0) {
              /* The point should be inside */
              SC_CHECK_ABORTF (point_is_inside,
                               "The point (%g,%g,%g) should be inside the unit quad, but isn't detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            else {
              /* The point should not be inside */
              SC_CHECK_ABORTF (!point_is_inside,
                               "The point (%g,%g,%g) should not be inside the unit quad, but is detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            break;
          case T8_ECLASS_TRIANGLE:
            /* The point is inside if all its x,y coordinates are in [0, 1], z = 0, and y <= x */
            if (test_point[0] >= 0 && test_point[1] >= 0 && test_point[2] == 0
                && test_point[0] <= 1 && test_point[1] <= test_point[0]) {
              /* The point should be inside */
              SC_CHECK_ABORTF (point_is_inside,
                               "The point (%g,%g,%g) should be inside the triangle, but isn't detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            else {
              /* The point should not be inside */
              SC_CHECK_ABORTF (!point_is_inside,
                               "The point (%g,%g,%g) should not be inside the triangle, but is detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            break;
          case T8_ECLASS_HEX:
            /* The point is inside if all its x,y,z coordinates are in [0, 1] */
            if (test_point[0] >= 0 && test_point[1] >= 0 && test_point[2] >= 0
                && test_point[0] <= 1 && test_point[1] <= 1
                && test_point[2] <= 1) {
              /* The point should be inside */
              SC_CHECK_ABORTF (point_is_inside,
                               "The point (%g,%g,%g) should be inside the unit hex, but isn't detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            else {
              /* The point should not be inside */
              SC_CHECK_ABORTF (!point_is_inside,
                               "The point (%g,%g,%g) should not be inside the unit hex, but is detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            break;
          case T8_ECLASS_TET:
            /* The point is inside if all its x,y, z coordinates are in [0, 1], and y <= z and x >= z */
            if (test_point[0] >= 0 && test_point[1] >= 0 && test_point[2] >= 0  /* x,y,z >= 0 */
                && test_point[0] <= 1   /* x <= 1 */
                && test_point[2] <= 1   /* z <= 1 */
                && test_point[1] <= test_point[2]       /* y <= z and x >= z */
                &&test_point[0] >= test_point[2]) {
              /* The point should be inside */
              SC_CHECK_ABORTF (point_is_inside,
                               "The point (%g,%g,%g) should be inside the tetrahedron, but isn't detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            else {
              /* The point should not be inside */
              SC_CHECK_ABORTF (!point_is_inside,
                               "The point (%g,%g,%g) should not be inside the tetrahedron, but is detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            break;
          case T8_ECLASS_PRISM:
            /* The point is inside if x and y are inside the triangle (x,y in [0,1], y <= x)
             * and if z in [0, 1] */
            if (test_point[0] >= 0 && test_point[1] >= 0 && test_point[2] >= 0
                && test_point[0] <= 1 && test_point[1] <= 1
                && test_point[2] <= 1 && test_point[1] <= test_point[0]) {
              /* The point should be inside */
              SC_CHECK_ABORTF (point_is_inside,
                               "The point (%g,%g,%g) should be inside the prism, but isn't detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            else {
              /* The point should not be inside */
              SC_CHECK_ABORTF (!point_is_inside,
                               "The point (%g,%g,%g) should not be inside the prism, but is detected.",
                               test_point[0], test_point[1], test_point[2]);
            }
            break;
          default:
            SC_ABORTF ("Point inside test not implemented for class %s.",
                       t8_eclass_to_string[eclass]);
          }
        }
      }
    }
  }

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

  /* TODO: Activate test for lines when point inside is implemented for lines */
  for (ieclass = T8_ECLASS_QUAD; ieclass < T8_ECLASS_COUNT; ieclass++) {
    if (ieclass != T8_ECLASS_PYRAMID) {
      /* TODO: does not work with pyramids yet */
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
