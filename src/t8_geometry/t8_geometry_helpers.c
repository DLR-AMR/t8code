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
t8_geom_linear_interpolation (const double *coefficients, const double *corner_values, const int corner_value_dim,
                              const int interpolation_dim, double *evaluated_function)
{
  double temp[3] = { 0 };
  for (int i_dim = 0; i_dim < corner_value_dim; i_dim++) {
    temp[i_dim] = corner_values[0 * corner_value_dim + i_dim] * (1 - coefficients[0]) /* x=0 y=0 z=0 */
                  + corner_values[1 * corner_value_dim + i_dim] * coefficients[0];    /* x=1 y=0 z=0 */
    if (interpolation_dim > 1) {
      temp[i_dim] *= (1 - coefficients[1]);
      temp[i_dim] += (corner_values[2 * corner_value_dim + i_dim] * (1 - coefficients[0]) /* x=0 y=1 z=0 */
                      + corner_values[3 * corner_value_dim + i_dim] * coefficients[0])    /* x=1 y=1 z=0 */
                     * coefficients[1];
      if (interpolation_dim == 3) {
        temp[i_dim] *= (1 - coefficients[2]);
        temp[i_dim]
          += (corner_values[4 * corner_value_dim + i_dim] * (1 - coefficients[0])
                * (1 - coefficients[1])                                                               /* x=0 y=0 z=1 */
              + corner_values[5 * corner_value_dim + i_dim] * coefficients[0] * (1 - coefficients[1]) /* x=1 y=0 z=1 */
              + corner_values[6 * corner_value_dim + i_dim] * (1 - coefficients[0]) * coefficients[1] /* x=0 y=1 z=1 */
              + corner_values[7 * corner_value_dim + i_dim] * coefficients[0] * coefficients[1])      /* x=1 y=1 z=1 */
             * coefficients[2];
      }
    }
    evaluated_function[i_dim] = temp[i_dim];
  }
}

void
t8_geom_triangular_interpolation (const double *coefficients, const double *corner_values, const int corner_value_dim,
                                  const int interpolation_dim, double *evaluated_function)
{
  /* The algorithm is able to calculate any point in a triangle or tetrahedron using barycentric coordinates.
   * All points are calculated by the sum of each corner point (e.g. p1 -> corner point 1) multiplied by a
   * scalar, which in this case are the reference coordinates (ref_coords).
   */
  double temp[3] = { 0 };

  for (int i_dim = 0; i_dim < corner_value_dim; i_dim++) {
    temp[i_dim] = (corner_values[corner_value_dim + i_dim] - /* (p2 - p1) * ref_coords */
                   corner_values[i_dim])
                    * coefficients[0]
                  + (interpolation_dim == 3
                       ? (corner_values[3 * corner_value_dim + i_dim] - corner_values[2 * corner_value_dim + i_dim])
                           * coefficients[1]
                       : 0.) /* (p4 - p3) * ref_coords */
                  + (corner_values[2 * corner_value_dim + i_dim] - corner_values[corner_value_dim + i_dim])
                      * coefficients[interpolation_dim - 1] /* (p3 - p2) * ref_coords */
                  + corner_values[i_dim];                   /* p1 */
    evaluated_function[i_dim] = temp[i_dim];
  }
}

void
t8_geom_compute_linear_geometry (t8_eclass_t tree_class, const double *tree_vertices, const double *ref_coords,
                                 const size_t num_coords, double *out_coords)
{
  int i_dim;
  size_t i_coord;
  const int dimension = t8_eclass_to_dimension[tree_class];
  /* Compute the coordinates, depending on the shape of the element */
  switch (tree_class) {
  case T8_ECLASS_VERTEX:
    /* A vertex has exactly one corner, and we already know its coordinates, since they are
     * the same as the trees coordinates. */
    for (i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
      for (i_dim = 0; i_dim < T8_ECLASS_MAX_DIM; i_dim++) {
        out_coords[offset_domain_dim + i_dim] = tree_vertices[offset_domain_dim + i_dim];
      }
    }
    break;
  case T8_ECLASS_TRIANGLE:
  case T8_ECLASS_TET:
    for (i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_tree_dim = i_coord * dimension;
      const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
      t8_geom_triangular_interpolation (ref_coords + offset_tree_dim, tree_vertices, T8_ECLASS_MAX_DIM, dimension,
                                        out_coords + offset_domain_dim);
    }
    break;
  case T8_ECLASS_PRISM: {
    double tri_vertices[9];
    double line_vertices[6];
    for (i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_tree_dim = i_coord * dimension;
      const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
      /* Prisminterpolation, via height and triangle */
      /* Get a triangle at the specific height */
      for (int i_tri_vertex = 0; i_tri_vertex < 3; i_tri_vertex++) {
        /* Vertices of each edge have to be linear in memory */
        memcpy (line_vertices, tree_vertices + i_tri_vertex * T8_ECLASS_MAX_DIM, T8_ECLASS_MAX_DIM * sizeof (double));
        memcpy (line_vertices + 3, tree_vertices + (i_tri_vertex + 3) * T8_ECLASS_MAX_DIM,
                T8_ECLASS_MAX_DIM * sizeof (double));
        t8_geom_linear_interpolation (ref_coords + offset_tree_dim + 2, line_vertices, T8_ECLASS_MAX_DIM, 1,
                                      tri_vertices + i_tri_vertex * T8_ECLASS_MAX_DIM);
      }
      t8_geom_triangular_interpolation (ref_coords + offset_tree_dim, tri_vertices, T8_ECLASS_MAX_DIM, 2,
                                        out_coords + offset_domain_dim);
    }
  } break;
  case T8_ECLASS_LINE:
  case T8_ECLASS_QUAD:
  case T8_ECLASS_HEX:
    for (i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_tree_dim = i_coord * dimension;
      const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
      t8_geom_linear_interpolation (ref_coords + offset_tree_dim, tree_vertices, T8_ECLASS_MAX_DIM, dimension,
                                    out_coords + offset_domain_dim);
    }
    break;
  case T8_ECLASS_PYRAMID: {
    double base_coords[2];
    double vec[3];
    for (i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_tree_dim = i_coord * dimension;
      const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
      /* Pyramid interpolation. After projecting the point onto the base, we use a bilinear interpolation to do a quad
       * interpolation on the base and then we interpolate via the height to the top vertex  */

      /* Project point on base */
      if (ref_coords[offset_tree_dim + 2] != 1.) {
        for (i_dim = 0; i_dim < 2; i_dim++) {
          base_coords[i_dim] = 1 - (1 - ref_coords[offset_tree_dim + i_dim]) / (1 - ref_coords[offset_tree_dim + 2]);
        }
      }
      else {
        for (i_dim = 0; i_dim < T8_ECLASS_MAX_DIM; i_dim++) {
          out_coords[offset_domain_dim + i_dim] = tree_vertices[4 * T8_ECLASS_MAX_DIM + i_dim];
        }
        continue;
      }
      /* Get a quad interpolation of the base */
      t8_geom_linear_interpolation (base_coords, tree_vertices, T8_ECLASS_MAX_DIM, 2, out_coords + offset_domain_dim);
      /* Get vector from base to pyramid tip */
      t8_vec_diff (tree_vertices + 4 * T8_ECLASS_MAX_DIM, out_coords + offset_domain_dim, vec);
      /* Add vector to base */
      for (i_dim = 0; i_dim < 3; i_dim++) {
        out_coords[offset_domain_dim + i_dim] += vec[i_dim] * ref_coords[offset_tree_dim + 2];
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
t8_geom_compute_linear_axis_aligned_geometry (t8_eclass_t tree_class, const double *tree_vertices,
                                              const double *ref_coords, const size_t num_coords, double out_coords[3])
{
  if (tree_class != T8_ECLASS_LINE && tree_class != T8_ECLASS_QUAD && tree_class != T8_ECLASS_HEX) {
    SC_ABORT ("Linear geometry coordinate computation is only supported for lines/quads/hexes.");
  }
#if T8_ENABLE_DEBUG
  /* Check if vertices are axis-aligned */
  if (tree_class == T8_ECLASS_LINE || tree_class == T8_ECLASS_QUAD) {
    /* The two vertices of a line must have two matching coordinates to be
     * axis-aligned. A quad needs one matching coordinate. */
    int n_equal_coords = 0;
    for (int i_dim = 0; i_dim < T8_ECLASS_MAX_DIM; ++i_dim) {
      if (abs (tree_vertices[i_dim] - tree_vertices[T8_ECLASS_MAX_DIM + i_dim]) <= SC_EPS) {
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
  const int dimension = t8_eclass_to_dimension[tree_class];
  /* Compute vector between both points */
  double vector[3];
  t8_vec_diff (tree_vertices + T8_ECLASS_MAX_DIM, tree_vertices, vector);

  /* Compute the coordinates of the reference point. */
  for (size_t i_coord = 0; i_coord < num_coords; ++i_coord) {
    const size_t offset_tree_dim = i_coord * dimension;
    const size_t offset_domain_dim = i_coord * T8_ECLASS_MAX_DIM;
    for (int i_dim = 0; i_dim < T8_ECLASS_MAX_DIM; ++i_dim) {
      out_coords[offset_domain_dim + i_dim] = tree_vertices[i_dim];
      out_coords[offset_domain_dim + i_dim] += ref_coords[offset_tree_dim] * vector[i_dim];
    }
  }
}

void
t8_geom_get_face_vertices (const t8_eclass_t tree_class, const double *tree_vertices, int face_index, int dim,
                           double *face_vertices)
{
  const int face_class = t8_eclass_face_types[tree_class][face_index];
  for (int i_face_vertex = 0; i_face_vertex < t8_eclass_num_vertices[face_class]; ++i_face_vertex) {
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
      const int i_tree_vertex = t8_face_vertex_to_tree_vertex[tree_class][face_index][i_face_vertex];
      face_vertices[i_face_vertex * dim + i_dim] = tree_vertices[i_tree_vertex * dim + i_dim];
    }
  }
}

void
t8_geom_get_edge_vertices (const t8_eclass_t tree_class, const double *tree_vertices, int edge_index, int dim,
                           double *edge_vertices)
{
  T8_ASSERT (t8_eclass_to_dimension[tree_class] == 3);
  for (int i_edge_vertex = 0; i_edge_vertex < 2; ++i_edge_vertex) {
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
      const int i_tree_vertex = t8_edge_vertex_to_tree_vertex[tree_class][edge_index][i_edge_vertex];
      edge_vertices[i_edge_vertex * dim + i_dim] = tree_vertices[i_tree_vertex * dim + i_dim];
    }
  }
}

void
t8_geom_get_ref_intersection (int edge_index, const double *ref_coords, double ref_intersection[2])
{
  double ref_slope;
  const t8_eclass_t eclass = T8_ECLASS_TRIANGLE;
  /* The opposite vertex of an edge always has the same index as the edge (see picture below). */
  const double *ref_opposite_vertex = t8_element_corner_ref_coords[eclass][edge_index];
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
    ref_slope = (ref_opposite_vertex[1] - ref_coords[1]) / (ref_opposite_vertex[0] - ref_coords[0]);
  }
  /* Now that we have the slope of the lines going through each vertex and the reference point,
   * we can calculate the intersection of the lines with each edge.
   */
  switch (edge_index) {
  case 0: /* edge 0 */
    /* Because of the verticality of edge 0, the x value of the intersection is always 1.
     * The y value is determined by the slope times the horizontal edge length.
     */
    ref_intersection[0] = 1;
    ref_intersection[1] = ref_slope * 1;
    break;
  case 1: /* edge 1 */
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
      ref_intersection[0] = ref_intersection[1]
        = ((ref_coords[0] * ref_opposite_vertex[1] - ref_coords[1] * ref_opposite_vertex[0])
             / -(ref_coords[1] - ref_opposite_vertex[1])
           + (ref_coords[0] - ref_opposite_vertex[0]));
      break;
    }
  case 2: /* edge 2 */
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
      ref_intersection[0] = (ref_coords[0] * ref_opposite_vertex[1] - ref_coords[1] * ref_opposite_vertex[0])
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
t8_geom_get_triangle_scaling_factor (int edge_index, const double *tree_vertices, const double *glob_intersection,
                                     const double *glob_ref_point)
{
  double dist_intersection, dist_ref;
  /* The scaling factor depends on the relation of the distance of the opposite vertex
   * to the global reference point and the distance of the opposite vertex to the global
   * intersection on the edge.
   */
  dist_intersection = sqrt (
    ((tree_vertices[edge_index * 3] - glob_intersection[0]) * (tree_vertices[edge_index * 3] - glob_intersection[0]))
    + ((tree_vertices[edge_index * 3 + 1] - glob_intersection[1])
       * (tree_vertices[edge_index * 3 + 1] - glob_intersection[1]))
    + ((tree_vertices[edge_index * 3 + 2] - glob_intersection[2])
       * (tree_vertices[edge_index * 3 + 2] - glob_intersection[2])));
  dist_ref
    = sqrt (((tree_vertices[edge_index * 3] - glob_ref_point[0]) * (tree_vertices[edge_index * 3] - glob_ref_point[0]))
            + ((tree_vertices[edge_index * 3 + 1] - glob_ref_point[1])
               * (tree_vertices[edge_index * 3 + 1] - glob_ref_point[1]))
            + ((tree_vertices[edge_index * 3 + 2] - glob_ref_point[2])
               * (tree_vertices[edge_index * 3 + 2] - glob_ref_point[2])));
  /* The closer the reference point is to the intersection, the bigger is the scaling factor. */
  double scaling_factor = dist_ref / dist_intersection;
  return scaling_factor;
}

double
t8_geom_get_scaling_factor_of_edge_on_face (const int edge, const int face, const double *ref_coords)
{
  /* Save the orthogonal direction and the maximum of that direction
   * of a tetrahedron edge in reference space on one of the neighbouring faces. 
   *           /|
   *          / |
   *         /  |
   *        /   |
   *       /    |--edge
   *      /   --|--face
   *     /      |
   *    /    <~~| orthogonal direction
   *   /<----o--| maximum othogonal direction
   *  /_________|*/
  double orthogonal_direction;
  double max_orthogonal_direction;
  double scaling_factor;

  switch (edge) { /* Check for edge of tetrahedron */
  case 0:
    switch (face) { /* Check for neighbouring face of edge 0 */
    case 2:
      orthogonal_direction = ref_coords[1];
      max_orthogonal_direction = ref_coords[0];
      break;
    case 3:
      orthogonal_direction = ref_coords[2];
      max_orthogonal_direction = ref_coords[0];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  case 1:
    switch (face) { /* Check for neighbouring face of edge 1 */
    case 1:
      orthogonal_direction = ref_coords[1];
      max_orthogonal_direction = ref_coords[0];
      break;
    case 3:
      orthogonal_direction = ref_coords[0] - ref_coords[2];
      max_orthogonal_direction = ref_coords[0];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  case 2:
    switch (face) { /* Check for neighbouring face of edge 2 */
    case 1:
      orthogonal_direction = ref_coords[0] - ref_coords[1];
      max_orthogonal_direction = ref_coords[0];
      break;
    case 2:
      orthogonal_direction = ref_coords[0] - ref_coords[2];
      max_orthogonal_direction = ref_coords[0];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  case 3:
    switch (face) { /* Check for neighbouring face of edge 3 */
    case 0:
      orthogonal_direction = ref_coords[1];
      max_orthogonal_direction = ref_coords[2];
      break;
    case 3:
      orthogonal_direction = 1 - ref_coords[0];
      max_orthogonal_direction = 1 - ref_coords[2];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  case 4:
    switch (face) { /* Check for neighbouring face of edge 4 */
    case 0:
      orthogonal_direction = ref_coords[2] - ref_coords[1];
      max_orthogonal_direction = ref_coords[2];
      break;
    case 2:
      orthogonal_direction = 1 - ref_coords[0];
      max_orthogonal_direction = 1 - ref_coords[2];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  case 5:
    switch (face) { /* Check for neighbouring face of edge 5 */
    case 0:
      orthogonal_direction = 1 - ref_coords[2];
      max_orthogonal_direction = 1 - ref_coords[1];
      break;
    case 1:
      orthogonal_direction = 1 - ref_coords[0];
      max_orthogonal_direction = 1 - ref_coords[1];
      break;
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }

  /* If the maximum orthogonal direction is 0 or 1, the reference coordinate lies on
   * one of the edge nodes and the scaling factor is therefore 0, because the displacement
   * at the nodes is always 0.
   * In all other cases the scaling factor is determined with one minus the relation of the orthogonal direction
   * to the maximum orthogonal direction. */
  if (max_orthogonal_direction == 0 || max_orthogonal_direction == 1) {
    scaling_factor = 0;
  }
  else {
    scaling_factor = 1 - (orthogonal_direction / max_orthogonal_direction);
  }
  return scaling_factor;
}

void
t8_geom_get_tet_face_intersection (const int face, const double *ref_coords, double face_intersection[3])
{
  /* Save reference corner coordinates of the current face */
  double ref_face_vertex_coords[9];
  for (int i_face_vertex = 0; i_face_vertex < 3; ++i_face_vertex) {
    for (int dim = 0; dim < 3; ++dim) {
      const int i_tree_vertex = t8_face_vertex_to_tree_vertex[T8_ECLASS_TET][face][i_face_vertex];
      ref_face_vertex_coords[i_face_vertex * 3 + dim] = t8_element_corner_ref_coords[T8_ECLASS_TET][i_tree_vertex][dim];
    }
  }

  /* Save the opposite vertex of the face in reference space.
   * Opposite vertex of a face has the same index as the face. */
  const double *ref_opposite_vertex = t8_element_corner_ref_coords[T8_ECLASS_TET][face];

  /* Save the normal of the current face */
  double normal[3];
  for (int dim = 0; dim < 3; ++dim) {
    normal[dim] = t8_face_normal_tet[face][dim];
  }

  /* Calculate the vector from the opposite vertex to the
   * reference coordinate in reference space */
  double vector[3];
  for (int dim = 0; dim < 3; ++dim) {
    vector[dim] = ref_coords[dim] - ref_opposite_vertex[dim];
  }

  /* Calculate t to get the point on the ray (extension of vector), which lies on the face.
   * The vector will later be multiplied by t to get the exact distance from the opposite vertex to the face intersection. 
   * t = ((point on face - point on vector) * normal of face) / (vector * normal of face) */
  double denominator = 0;
  double numerator = 0;
  for (int dim = 0; dim < 3; ++dim) {
    denominator += (ref_face_vertex_coords[dim] - ref_opposite_vertex[dim]) * normal[dim];
    numerator += vector[dim] * normal[dim];
  }
  double t = denominator / numerator;

  /* Calculate face intersection by scaling vector with t.
   * If the reference coordinate is equal to the opposite vertex,
   * the intersection is equal to one of the ref_face_vertex_coords. */
  if (ref_coords[0] == ref_opposite_vertex[0] && ref_coords[1] == ref_opposite_vertex[1]
      && ref_coords[2] == ref_opposite_vertex[2]) {
    for (int dim = 0; dim < 3; ++dim) {
      face_intersection[dim] = ref_face_vertex_coords[dim];
    }
  }
  else {
    for (int dim = 0; dim < 3; ++dim) {
      face_intersection[dim] = ref_opposite_vertex[dim] + vector[dim] * t;
    }
  }
}
