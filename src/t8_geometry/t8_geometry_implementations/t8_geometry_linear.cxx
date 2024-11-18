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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_vec.h>

t8_geometry_linear::t8_geometry_linear (): t8_geometry_with_vertices ("t8_geom_linear")
{
}

t8_geometry_linear::~t8_geometry_linear ()
{
}

void
t8_geometry_linear::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                      const size_t num_coords, double *out_coords) const
{
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);
}

void
t8_geometry_linear::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                               const size_t num_coords, double *jacobian) const
{
  SC_ABORT ("Not implemented.");
}

#if T8_ENABLE_DEBUG
/* Test whether four given points in 3D are coplanar up to a given tolerance.
 */
static int
t8_four_points_coplanar (const double p_0[3], const double p_1[3], const double p_2[3], const double p_3[3],
                         const double tolerance)
{
  /* Let p0, p1, p2, p3 be the four points.
   * The four points are coplanar if the normal vectors to the triangles
   * p0, p1, p2 and p0, p2, p3 are pointing in the same direction.
   *
   * We build the vectors A = p1 - p0, B = p2 - p0 and C = p3 - p0.
   * The normal vectors to the triangles are n1 = A x B and n2 = A x C.
   * These are pointing in the same direction if their cross product is 0.
   * Hence we check if || n1 x n2 || < tolerance. */

  /* A = p1 - p0 */
  double A[3];
  t8_vec_axpyz (p_0, p_1, A, -1);

  /* B = p2 - p0 */
  double B[3];
  t8_vec_axpyz (p_0, p_2, B, -1);

  /* C = p3 - p0 */
  double C[3];
  t8_vec_axpyz (p_0, p_3, C, -1);

  /* n1 = A x B */
  double A_cross_B[3];
  t8_vec_cross (A, B, A_cross_B);

  /* n2 = A x C */
  double A_cross_C[3];
  t8_vec_cross (A, C, A_cross_C);

  /* n1 x n2 */
  double n1_cross_n2[3];
  t8_vec_cross (A_cross_B, A_cross_C, n1_cross_n2);

  /* || n1 x n2 || */
  const double norm = t8_vec_norm (n1_cross_n2);
  return norm < tolerance;
}
#endif

void
t8_geometry_linear::t8_geom_point_batch_inside_element (t8_forest_t forest, t8_locidx_t ltreeid,
                                                        const t8_element_t *element, const double *points,
                                                        const int num_points, int *is_inside,
                                                        const double tolerance) const
{
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_shape_t element_shape = ts->t8_element_shape (element);
  switch (element_shape) {
  case T8_ECLASS_VERTEX: {
    /* A point is 'inside' a vertex if they have the same coordinates */
    double vertex_coords[3];
    /* Get the vertex coordinates */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, vertex_coords);
    /* Check whether the point and the vertex are within tolerance distance
       * to each other */
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      is_inside[ipoint] = t8_vertex_point_inside (vertex_coords, &points[ipoint * 3], tolerance);
    }
    return;
  }
  case T8_ECLASS_LINE: {
    /* A point p is inside a line that is defined by the edge nodes
     * p_0 and p_1
     * if and only if the linear system
     * (p_1 - p_0)x = p - p_0
     * has a solution x with 0 <= x <= 1
     */
    double p_0[3], v[3];

    /* Compute the vertex coordinates of the line */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, p_0);
    /* v = p_1 */
    t8_forest_element_coordinate (forest, ltreeid, element, 1, v);
    /* v = p_1 - p_0 */
    t8_vec_axpy (p_0, v, -1);
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      is_inside[ipoint] = t8_line_point_inside (p_0, v, &points[ipoint * 3], tolerance);
    }
    return;
  }
  case T8_ECLASS_QUAD: {
    /* We divide the quad in two triangles and use the triangle check. */
    double p_0[3], p_1[3], p_2[3], p_3[3];
    /* Compute the vertex coordinates of the quad */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, p_0);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, p_1);
    t8_forest_element_coordinate (forest, ltreeid, element, 2, p_2);
    t8_forest_element_coordinate (forest, ltreeid, element, 3, p_3);

#if T8_ENABLE_DEBUG
    /* Issue a warning if the points of the quad do not lie in the same plane */
    if (!t8_four_points_coplanar (p_0, p_1, p_2, p_3, tolerance)) {
      t8_debugf ("WARNING: Testing if point is inside a quad that is not coplanar. This test will be inaccurate.\n");
    }
#endif
    double v[3];
    double w[3];
    /* v = v - p_0 = p_1 - p_0 */
    t8_vec_axpyz (p_0, p_1, v, -1);
    /* w = w - p_0 = p_2 - p_0 */
    t8_vec_axpyz (p_0, p_2, w, -1);
    /* Check whether the point is inside the first triangle. */
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      is_inside[ipoint] = t8_triangle_point_inside (p_0, v, w, &points[ipoint * 3], tolerance);
    }
    /* If not, check whether the point is inside the second triangle. */
    /* v = v - p_0 = p_1 - p_0 */
    t8_vec_axpyz (p_1, p_2, v, -1);
    /* w = w - p_0 = p_2 - p_0 */
    t8_vec_axpyz (p_1, p_3, w, -1);
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      if (!is_inside[ipoint]) {
        /* point_inside is true if the point was inside the first or second triangle. Otherwise it is false. */
        is_inside[ipoint] = t8_triangle_point_inside (p_1, v, w, &points[ipoint * 3], tolerance);
      }
    }
    return;
  }
  case T8_ECLASS_TRIANGLE: {
    double p_0[3], p_1[3], p_2[3];

    /* Compute the vertex coordinates of the triangle */
    t8_forest_element_coordinate (forest, ltreeid, element, 0, p_0);
    t8_forest_element_coordinate (forest, ltreeid, element, 1, p_1);
    t8_forest_element_coordinate (forest, ltreeid, element, 2, p_2);
    double v[3];
    double w[3];
    /* v = v - p_0 = p_1 - p_0 */
    t8_vec_axpyz (p_0, p_1, v, -1);
    /* w = w - p_0 = p_2 - p_0 */
    t8_vec_axpyz (p_0, p_2, w, -1);

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      is_inside[ipoint] = t8_triangle_point_inside (p_0, v, w, &points[ipoint * 3], tolerance);
    }
    return;
  }
  case T8_ECLASS_TET:
  case T8_ECLASS_HEX:
  case T8_ECLASS_PRISM:
  case T8_ECLASS_PYRAMID: {
    /* For bilinearly interpolated volume elements, a point is inside an element
     * if and only if it lies on the inner side of each face.
     * The inner side is defined as the side where the outside normal vector does not
     * point to.
     * The point is on this inner side if and only if the scalar product of
     * a point on the plane minus the point with the outer normal of the face
     * is >= 0.
     *
     * In other words, let p be the point to check, n the outer normal and x a point
     * on the plane, then p is on the inner side if and only if
     *  <x - p, n> >= 0
     */

    const int num_faces = ts->t8_element_num_faces (element);
    /* Assume that every point is inside of the element */
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      is_inside[ipoint] = 1;
    }
    std::array<std::array<double, 3>, T8_ECLASS_MAX_FACES> face_normals;
    for (int iface = 0; iface < num_faces; ++iface) {
      /* Compute the outer normal n of the face */
      t8_forest_element_face_normal (forest, ltreeid, element, iface, face_normals[iface].data ());
    }
    std::array<std::array<double, 3>, T8_ECLASS_MAX_FACES> face_points;
    for (int iface = 0; iface < num_faces; ++iface) {
      const int afacecorner = ts->t8_element_get_face_corner (element, iface, 0);
      t8_forest_element_coordinate (forest, ltreeid, element, afacecorner, face_points[iface].data ());
    }
    for (int ipoint = 0; ipoint < num_points; ipoint++) {
      for (int iface = 0; iface < num_faces; ++iface) {
        if (t8_plane_point_inside (face_points[iface].data (), face_normals[iface].data (), &points[ipoint * 3]) == 0) {
          is_inside[ipoint] = 0;
          break;
        }
      }
    }
    return;
  }
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear.h.
 * Create a new geometry. */
t8_geometry_c *
t8_geometry_linear_new ()
{
  t8_geometry_linear *geom = new t8_geometry_linear ();
  return (t8_geometry_c *) geom;
}

void
t8_geometry_linear_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_LINEAR);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();
