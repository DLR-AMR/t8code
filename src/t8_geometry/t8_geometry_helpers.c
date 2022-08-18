/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>

void
t8_geom_linear_interpolation (const double *coefficients,
                              const double *corner_values,
                              int corner_value_dim,
                              int interpolation_dim,
                              double *evaluated_function)
{
  double              temp[3] = { 0 };
  for (int i = 0; i < corner_value_dim; i++) {
    temp[i] = corner_values[0 * corner_value_dim + i] * (1 - coefficients[0])   /* x=0 y=0 z=0 */
      +corner_values[1 * corner_value_dim + i] * coefficients[0];       /* x=1 y=0 z=0 */
    if (interpolation_dim > 1) {
      temp[i] *= (1 - coefficients[1]);
      temp[i] += (corner_values[2 * corner_value_dim + i] * (1 - coefficients[0])       /* x=0 y=1 z=0 */
                  +corner_values[3 * corner_value_dim + i] * coefficients[0])   /* x=1 y=1 z=0 */
        *coefficients[1];
      if (interpolation_dim == 3) {
        temp[i] *= (1 - coefficients[2]);
        temp[i] += (corner_values[4 * corner_value_dim + i] * (1 - coefficients[0]) * (1 - coefficients[1])     /* x=0 y=0 z=1 */
                    +corner_values[5 * corner_value_dim + i] * coefficients[0] * (1 - coefficients[1])  /* x=1 y=0 z=1 */
                    +corner_values[6 * corner_value_dim + i] * (1 - coefficients[0]) * coefficients[1]  /* x=0 y=1 z=1 */
                    +corner_values[7 * corner_value_dim + i] * coefficients[0] * coefficients[1])       /* x=1 y=1 z=1 */
          *coefficients[2];
      }
    }
    evaluated_function[i] = temp[i];
  }
}

void
t8_geom_compute_linear_geometry (t8_eclass_t tree_class,
                                 const double *tree_vertices,
                                 const double *ref_coords,
                                 double out_coords[3])
{
  int                 i;
  int                 dimension = t8_eclass_to_dimension[tree_class];
  /* Compute the coordinates, depending on the shape of the element */
  switch (tree_class) {
  case T8_ECLASS_VERTEX:
    /* A vertex has exactly one corner, and we already know its coordinates, since they are
     * the same as the trees coordinates. */
    for (i = 0; i < 3; i++) {
      out_coords[i] = tree_vertices[i];
    }
    break;
  case T8_ECLASS_LINE:
    for (i = 0; i < 3; i++) {
      out_coords[i] =
        (tree_vertices[3 + i] -
         tree_vertices[i]) * ref_coords[0] + tree_vertices[i];
    }
    break;
  case T8_ECLASS_TRIANGLE:
  case T8_ECLASS_TET:
    for (i = 0; i < 3; i++) {
      out_coords[i] =
        (tree_vertices[3 + i] -
         tree_vertices[i]) * ref_coords[0] + (dimension ==
                                              3
                                              ? (tree_vertices
                                                 [9 + i] -
                                                 tree_vertices[6 + i])
                                              * ref_coords[1]
                                              : 0.)
        + (tree_vertices[6 + i] -
           tree_vertices[3 + i]) * ref_coords[dimension - 1]
        + tree_vertices[i];
    }
    break;
  case T8_ECLASS_PRISM:
    {
      /* Prisminterpolation, via height, and triangle */
      /* Get a triangle at the specific height */
      double              tri_vertices[9];
      for (i = 0; i < 9; i++) {
        tri_vertices[i] =
          (tree_vertices[9 + i] -
           tree_vertices[i]) * ref_coords[2] + tree_vertices[i];
      }
      for (i = 0; i < 3; i++) {
        out_coords[i] =
          (tri_vertices[3 + i] - tri_vertices[i]) * ref_coords[0] +
          (tri_vertices[6 + i] - tri_vertices[3 + i]) * ref_coords[1]
          + tri_vertices[i];
      }
    }
    break;
  case T8_ECLASS_QUAD:
  case T8_ECLASS_HEX:
    t8_geom_linear_interpolation (ref_coords,
                                  tree_vertices, 3, dimension, out_coords);
    break;
  case T8_ECLASS_PYRAMID:
    {
      double              ray[3], lambda, quad_coords[3], length, length2;
      length = 0;
      length2 = 0;
      quad_coords[2] = 0;
      /*vertices = tree_vertices */
      /*coordinates = out_coords */
      /*vertex_coords = len * corner_coords = ref_coords */

      /*In this case, the vertex is the tip of the parent pyramid and we don't have to compute
       * anything.*/
      if (ref_coords[0] == 1. && ref_coords[1] == 1. && ref_coords[2] == 1.) {
        for (i = 0; i < 3; i++) {
          out_coords[i] = tree_vertices[12 + i];
        }
        break;
      }
      /* Project vertex_coord onto x-y-plane */
      for (i = 0; i < 3; i++) {
        ray[i] = 1 - ref_coords[i];
      }
      lambda = ref_coords[2] / ray[2];
      for (i = 0; i < 2; i++) {
        /*Compute coords of vertex in the plane */
        quad_coords[i] = ref_coords[i] - lambda * ray[i];
        length += (1 - quad_coords[i]) * (1 - quad_coords[i]);
      }
      length += 1;
      /*compute the ratio */
      for (i = 0; i < 3; i++) {
        length2 +=
          (ref_coords[i] - quad_coords[i]) * (ref_coords[i] - quad_coords[i]);
      }
      lambda = sqrt (length2) / sqrt (length);

      /*Interpolate on quad */
      t8_geom_linear_interpolation ((const double *) quad_coords,
                                    tree_vertices, 3, 2, out_coords);
      /*Project it back */
      for (i = 0; i < 3; i++) {
        out_coords[i] += (tree_vertices[12 + i] - out_coords[i]) * lambda;
      }
      break;
    }
  default:
    SC_ABORT ("Linear geometry coordinate computation is supported only for "
              "vertices/lines/triangles/tets/quads/prisms/hexes/pyramids.");
  }
}

void
t8_geom_get_face_vertices (const t8_eclass_t tree_class,
                           const double *tree_vertices,
                           int face_index, int dim, double *face_vertices)
{
  const int           face_class =
    t8_eclass_face_types[tree_class][face_index];
  for (int i_face_vertex = 0;
       i_face_vertex < t8_eclass_num_vertices[face_class]; ++i_face_vertex) {
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
      const int           i_tree_vertex =
        t8_face_vertex_to_tree_vertex[tree_class][face_index][i_face_vertex];
      face_vertices[i_face_vertex * dim + i_dim] =
        tree_vertices[i_tree_vertex * dim + i_dim];
    }
  }
}

void
t8_geom_get_edge_vertices (const t8_eclass_t tree_class,
                           const double *tree_vertices,
                           int edge_index, int dim, double *edge_vertices)
{
  T8_ASSERT (t8_eclass_to_dimension[tree_class] == 3);
  for (int i_edge_vertex = 0; i_edge_vertex < 2; ++i_edge_vertex) {
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
      const int           i_tree_vertex =
        t8_edge_vertex_to_tree_vertex[edge_index][i_edge_vertex];
      edge_vertices[i_edge_vertex * dim + i_dim] =
        tree_vertices[i_tree_vertex * dim + i_dim];
    }
  }
}
