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

#include <t8_vec.h>
#include <t8_eclass.h>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.h>

void
t8_geom_linear_interpolation (const double *coefficients,
                              const double *corner_values,
                              const int corner_value_dim,
                              const int interpolation_dim,
                              double *evaluated_function)
{
  double              temp[3] = { 0 };
  for (int dim = 0; dim < corner_value_dim; dim++) {
    temp[dim] = corner_values[0 * corner_value_dim + dim] * (1 - coefficients[0])       /* x=0 y=0 z=0 */
      +corner_values[1 * corner_value_dim + dim] * coefficients[0];     /* x=1 y=0 z=0 */
    if (interpolation_dim > 1) {
      temp[dim] *= (1 - coefficients[1]);
      temp[dim] += (corner_values[2 * corner_value_dim + dim] * (1 - coefficients[0])   /* x=0 y=1 z=0 */
                    +corner_values[3 * corner_value_dim + dim] * coefficients[0])       /* x=1 y=1 z=0 */
        *coefficients[1];
      if (interpolation_dim == 3) {
        temp[dim] *= (1 - coefficients[2]);
        temp[dim] += (corner_values[4 * corner_value_dim + dim] * (1 - coefficients[0]) * (1 - coefficients[1]) /* x=0 y=0 z=1 */
                      +corner_values[5 * corner_value_dim + dim] * coefficients[0] * (1 - coefficients[1])      /* x=1 y=0 z=1 */
                      +corner_values[6 * corner_value_dim + dim] * (1 - coefficients[0]) * coefficients[1]      /* x=0 y=1 z=1 */
                      +corner_values[7 * corner_value_dim + dim] * coefficients[0] * coefficients[1])   /* x=1 y=1 z=1 */
          *coefficients[2];
      }
    }
    evaluated_function[dim] = temp[dim];
  }
}

void
t8_geom_triangular_interpolation (const double *coefficients,
                                  const double *corner_values,
                                  const int corner_value_dim,
                                  const int interpolation_dim,
                                  double *evaluated_function)
{
  /* The algorithm is able to calculate any point in a triangle or tetrahedron using barycentric coordinates.
   * All points are calculated by the sum of each corner point (e.g. p1 -> corner point 1) multiplied by a
   * scalar, which in this case are the reference coordinates (ref_coords).
   */
  double              temp[3] = { 0 };

  for (int dim = 0; dim < corner_value_dim; dim++) {
    temp[dim] = (corner_values[corner_value_dim + dim] -        /* (p2 - p1) * ref_coords */
                 corner_values[dim]) * coefficients[0] + (interpolation_dim ==
                                                          3
                                                          ? (corner_values
                                                             [3 *
                                                              corner_value_dim
                                                              + dim] -
                                                             corner_values[2 *
                                                                           corner_value_dim
                                                                           +
                                                                           dim])
                                                          * coefficients[1]
                                                          : 0.) /* (p4 - p3) * ref_coords */
      +(corner_values[2 * corner_value_dim + dim] - corner_values[corner_value_dim + dim]) * coefficients[interpolation_dim - 1]        /* (p3 - p2) * ref_coords */
      +corner_values[dim];      /* p1 */
    evaluated_function[dim] = temp[dim];
  }
}

void
t8_geom_compute_linear_geometry (t8_eclass_t tree_class,
                                 const double *tree_vertices,
                                 const double *ref_coords,
                                 const int num_coords, double out_coords[3])
{
  int                 dim;
  int                 coord;
  const int           dimension = t8_eclass_to_dimension[tree_class];
  /* Compute the coordinates, depending on the shape of the element */
  switch (tree_class) {
  case T8_ECLASS_VERTEX:
    /* A vertex has exactly one corner, and we already know its coordinates, since they are
     * the same as the trees coordinates. */
    for (coord = 0; coord < num_coords; coord++) {
      const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
      for (dim = 0; dim < T8_ECLASS_MAX_DIM; dim++) {
        out_coords[offset_domain_dim + dim] =
          tree_vertices[offset_domain_dim + dim];
      }
    }
    break;
  case T8_ECLASS_TRIANGLE:
  case T8_ECLASS_TET:
    for (coord = 0; coord < num_coords; coord++) {
      const int           offset_tree_dim = coord * dimension;
      const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
      t8_geom_triangular_interpolation (ref_coords + offset_tree_dim,
                                        tree_vertices, T8_ECLASS_MAX_DIM,
                                        dimension,
                                        out_coords + offset_domain_dim);
    }
    break;
  case T8_ECLASS_PRISM:
    {
      for (coord = 0; coord < num_coords; coord++) {
        const int           offset_tree_dim = coord * dimension;
        const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
        /* Prisminterpolation, via height, and triangle */
        /* Get a triangle at the specific height */
        double              tri_vertices[9];
        for (dim = 0; dim < 3; dim++) {
          t8_geom_linear_interpolation (ref_coords + offset_tree_dim + 2,
                                        tree_vertices, T8_ECLASS_MAX_DIM, 1,
                                        tri_vertices);
        }
        t8_geom_triangular_interpolation (ref_coords + offset_tree_dim,
                                          tri_vertices, T8_ECLASS_MAX_DIM, 2,
                                          out_coords + offset_domain_dim);
      }
    }
    break;
  case T8_ECLASS_LINE:
  case T8_ECLASS_QUAD:
  case T8_ECLASS_HEX:
    for (coord = 0; coord < num_coords; coord++) {
      const int           offset_tree_dim = coord * dimension;
      const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
      t8_geom_linear_interpolation (ref_coords + offset_tree_dim,
                                    tree_vertices, T8_ECLASS_MAX_DIM,
                                    dimension,
                                    out_coords + offset_domain_dim);
    }
    break;
  case T8_ECLASS_PYRAMID:
    {
      for (coord = 0; coord < num_coords; coord++) {
        const int           offset_tree_dim = coord * dimension;
        const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
        double              base_coords[2];
        double              vec[3];
        /* Pyramid interpolation. After projecting the point onto the base,
         * we use a bilinear interpolation, to do a quad interpolation on the base
         * and then we interpolate via the height to the top vertex */
        /* Project point on base */
        for (dim = 0; dim < 2; dim++) {
          base_coords[dim] =
            ref_coords[offset_tree_dim + dim] / (1. -
                                                 ref_coords[offset_tree_dim +
                                                            2]);
        }
        /* Get a quad interpolation of the base */
        t8_geom_linear_interpolation (base_coords,
                                      tree_vertices, T8_ECLASS_MAX_DIM, 2,
                                      out_coords + offset_domain_dim);
        /* Get vector from base to pyramid tip */
        t8_vec_diff (out_coords + offset_domain_dim,
                     tree_vertices + 4 * T8_ECLASS_MAX_DIM, vec);
        /* Linear interpolation on vector. */
        t8_geom_linear_interpolation (ref_coords + offset_tree_dim + 2, vec,
                                      T8_ECLASS_MAX_DIM, 1, vec);
        /* Add vector to base */
        for (dim = 0; dim < 3; dim++) {
          out_coords[offset_domain_dim + dim] += vec[dim];
        }
      }
      break;
    }
  default:
    SC_ABORT ("Linear geometry coordinate computation is only supported for "
              "vertices/lines/triangles/tets/quads/prisms/hexes/pyramids.");
  }
}

void
t8_geom_compute_linear_axis_aligned_geometry (t8_eclass_t tree_class,
                                              const double *tree_vertices,
                                              const double *ref_coords,
                                              const int num_coords,
                                              double out_coords[3])
{
  if (tree_class != T8_ECLASS_LINE && tree_class != T8_ECLASS_QUAD
      && tree_class != T8_ECLASS_HEX) {
    SC_ABORT ("Linear geometry coordinate computation is only supported for "
              "lines/quads/hexes.");
  }
#if T8_ENABLE_DEBUG
  /* Check if vertices are axis-aligned */
  if (tree_class == T8_ECLASS_LINE || tree_class == T8_ECLASS_QUAD) {
    /* The two vertices of a line must have two matching coordinates to be
     * axis-aligned. A quad needs one matching coordinate. */
    int                 n_equal_coords = 0;
    for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
      if (abs (tree_vertices[dim] - tree_vertices[T8_ECLASS_MAX_DIM + dim]) <=
          SC_EPS) {
        ++n_equal_coords;
      }
    }
    if (tree_class == T8_ECLASS_LINE && n_equal_coords != 2) {
      SC_ABORT ("Line vertices are not axis-aligned.");
    }
    else if (tree_class == T8_ECLASS_QUAD && n_equal_coords != 1) {
      SC_ABORT ("Quad vertices are not axis-aligned.");
    }
  }
#endif /* T8_ENABLE_DEBUG */
  const int           dimension = t8_eclass_get_dimension (tree_class);
  /* Compute vector between both points */
  double              vector[3];
  t8_vec_diff (tree_vertices + T8_ECLASS_MAX_DIM, tree_vertices, vector);

  /* Compute the coordinates of the reference point. */
  for (int coord = 0; coord < num_coords; ++coord) {
    const int           offset_tree_dim = coord * dimension;
    const int           offset_domain_dim = coord * T8_ECLASS_MAX_DIM;
    for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
      out_coords[offset_tree_dim + dim] = tree_vertices[dim];
      out_coords[offset_tree_dim + dim] +=
        ref_coords[offset_tree_dim] * vector[dim];
    }
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

void
t8_geom_get_ref_intersection (int edge_index,
                              const double *ref_coords,
                              double ref_intersection[2])
{
  double              ref_slope;
  const t8_eclass_t   eclass = T8_ECLASS_TRIANGLE;
  /* The opposite vertex of an edge always has the same index as the edge (see picture below). */
  const double       *ref_opposite_vertex =
    t8_element_corner_ref_coords[eclass][edge_index];
  /*              2
   *            / |
   *           /  |
   *          /   |
   *         /    |
   *        E1   E0
   *       /      |
   *      /       |
   *     /        |
   *    /         |
   *   0----E2----1
   *
   * First, we calculate the slope of the line going through the reference point
   * and the opposite vertex for each edge of the triangle.
   */

  /* In case the reference point is equal to the opposite vertex, the slope of the line is 0. */
  if (ref_opposite_vertex[0] == ref_coords[0]) {
    ref_slope = 0;
  }
  /* slope = (y2-y1)/(x2-x1) */
  else {
    ref_slope = (ref_opposite_vertex[1] - ref_coords[1])
      / (ref_opposite_vertex[0] - ref_coords[0]);
  }
  /* Now that we have the slope of the lines going through each vertex and the reference point,
   * we can calculate the intersection of the lines with each edge.
   */
  switch (edge_index) {
  case 0:                      /* edge 0 */
    /* Because of the verticality of edge 0, the x value of the intersection is always 1.
     * The y value is determined by the slope times the horizontal edge length.
     */
    ref_intersection[0] = 1;
    ref_intersection[1] = ref_slope * 1;
    break;
  case 1:                      /* edge 1 */
    /* If the reference point lies somewhere on edge 0, the intersection has to be at (1,1). */
    if (ref_coords[0] == ref_opposite_vertex[0]) {
      ref_intersection[0] = 1;
      ref_intersection[1] = 1;
      break;
    }
    /* If the reference point lies somewhere on edge 2, the intersection has to be at (0,0). */
    else if (ref_coords[1] == ref_opposite_vertex[1]) {
      ref_intersection[0] = 0;
      ref_intersection[1] = 0;
      break;
    }
    else {
      /* intersectionX = (x1y2-y1x2)(x3-x4)-(x1-x2)(x3y4-y3x4)
       *                 /(x1-x2)(y3-y4)-(y1-y2)(x3-x4)
       * intersectionY = (x1y2-y1x2)(y3-y4)-(y1-y2)(x3y4-y3x4)
       *                 /(x1-x2)(y3-y4)-(y1-y2)(x3-x4)
       * 
       * x1=0 y1=0 x2=1 y2=1 x3=ref_coords[0] y3=ref_coords[1] x4=ref_opposite_vertex[0] y4=ref_opposite_vertex[1]
       * 
       * Since the intersection point lies on edge 2, which has a slope of 1, the x and the y value has to be equal
       */
      ref_intersection[0] = ref_intersection[1] =
        ((ref_coords[0] * ref_opposite_vertex[1] -
          ref_coords[1] * ref_opposite_vertex[0])
         / -(ref_coords[1] - ref_opposite_vertex[1]) +
         (ref_coords[0] - ref_opposite_vertex[0]));
      break;
    }
  case 2:                      /* edge 2 */
    /* If the reference point lies somewhere on edge 0, the intersection has to be at (1,0). */
    if (ref_coords[0] == ref_opposite_vertex[0]) {
      ref_intersection[0] = 1;
      ref_intersection[1] = 0;
      break;
    }
    /* If the reference point is equal to the opposite vertex, the intersection has to be at (0,1). */
    else if (ref_coords[1] == ref_opposite_vertex[1]) {
      ref_intersection[0] = 0;
      ref_intersection[1] = 1;
    }
    else {
      /* intersectionX = (x1y2-y1x2)(x3-x4)-(x1-x2)(x3y4-y3x4)
       *                 /(x1-x2)(y3-y4)-(y1-y2)(x3-x4)
       * intersectionY = (x1y2-y1x2)(y3-y4)-(y1-y2)(x3y4-y3x4)
       *                 /(x1-x2)(y3-y4)-(y1-y2)(x3-x4)
       * 
       * x1=0 y1=0 x2=1 y2=0 x3=ref_coords[0] y3=ref_coords[1] x4=ref_opposite_vertex[0] y4=ref_opposite_vertex[1]
       * 
       * Since the intersection point lies on edge 2, which has a slope of 1 in the reference space, the x and the y value has to be equal
       */
      ref_intersection[0] =
        (ref_coords[0] * ref_opposite_vertex[1] -
         ref_coords[1] * ref_opposite_vertex[0])
        / (-(ref_coords[1] - ref_opposite_vertex[1]));
      /* Since edge 2 is horizontal, the y value of the intersection always has to be 0. */
      ref_intersection[1] = 0;
      break;
    }
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
}

double
t8_geom_get_triangle_scaling_factor (int edge_index,
                                     const double *tree_vertices,
                                     const double *glob_intersection,
                                     const double *glob_ref_point)
{
  double              dist_intersection, dist_ref;
  /* The scaling factor depends on the relation of the distance of the opposite vertex
   * to the global reference point and the distance of the opposite vertex to the global
   * intersection on the edge.
   */
  dist_intersection =
    sqrt (((tree_vertices[edge_index * 3] -
            glob_intersection[0]) * (tree_vertices[edge_index * 3] -
                                     glob_intersection[0]))
          +
          ((tree_vertices[edge_index * 3 + 1] -
            glob_intersection[1]) * (tree_vertices[edge_index * 3 + 1] -
                                     glob_intersection[1]))
          +
          ((tree_vertices[edge_index * 3 + 2] -
            glob_intersection[2]) * (tree_vertices[edge_index * 3 + 2] -
                                     glob_intersection[2])));
  dist_ref =
    sqrt (((tree_vertices[edge_index * 3] -
            glob_ref_point[0]) * (tree_vertices[edge_index * 3] -
                                  glob_ref_point[0]))
          +
          ((tree_vertices[edge_index * 3 + 1] -
            glob_ref_point[1]) * (tree_vertices[edge_index * 3 + 1] -
                                  glob_ref_point[1]))
          +
          ((tree_vertices[edge_index * 3 + 2] -
            glob_ref_point[2]) * (tree_vertices[edge_index * 3 + 2] -
                                  glob_ref_point[2])));
  /* The closer the reference point is to the intersection, the bigger is the scaling factor. */
  double              scaling_factor = dist_ref / dist_intersection;
  return scaling_factor;
}
