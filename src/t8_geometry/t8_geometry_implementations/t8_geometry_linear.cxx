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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>

t8_geometry_linear::t8_geometry_linear (int dim): t8_geometry_with_vertices (dim, "")
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t num_chars = 100;
  char *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_linear_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_linear::~t8_geometry_linear ()
{
  T8_FREE ((char *) name);
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

/** Check if a point lies inside a vertex
 * 
 * \param[in] vertex_coords The coordinates of the vertex
 * \param[in] point         The coordinates of the point to check
 * \param[in] tolerance     A double > 0 defining the tolerance
 * \return                  0 if the point is outside, 1 otherwise.  
 */
static int
t8_vertex_point_inside (const double vertex_coords[3], const double point[3], const double tolerance)
{
  T8_ASSERT (tolerance > 0);
  if (t8_vec_dist (vertex_coords, point) > tolerance) {
    return 0;
  }
  return 1;
}

/**
 * Check if a point is inside a line that is defined by a starting point \a p_0
 * and a vector \a vec
 * 
 * \param[in] p_0         Starting point of the line
 * \param[in] vec         Direction of the line (not normalized)
 * \param[in] point       The coordinates of the point to check
 * \param[in] tolerance   A double > 0 defining the tolerance
 * \return                0 if the point is outside, 1 otherwise.  
 */
static int
t8_line_point_inside (const double *p_0, const double *vec, const double *point, const double tolerance)
{
  T8_ASSERT (tolerance > 0);
  double b[3];
  /* b = p - p_0 */
  t8_vec_axpyz (p_0, point, b, -1);
  double x = 0; /* Initialized to prevent compiler warning. */
  int i;
  /* So x is the solution to
  * vec * x = b.
  * We can compute it as
  * x = b[i] / vec[i]
  * if any vec[i] is not 0.
  *
  * Otherwise the line is degenerated (which should not happen).
  */
  for (i = 0; i < 3; ++i) {
    if (vec[i] != 0) {
      x = b[i] / vec[i];
      break; /* found a non-zero coordinate. We can stop now. */
    }
  }

  /* If i == 3 here, then vec = 0 and hence the line is degenerated. */
  SC_CHECK_ABORT (i < 3, "Degenerated line element. Both endpoints are the same.");

  if (x < -tolerance || x > 1 + tolerance) {
    /* x is not an admissible solution. */
    return 0;
  }

  /* we can check whether x gives us a solution by
     * checking whether
     *  vec * x = b
     * is actually true.
     */
  double vec_check[3] = { vec[0], vec[1], vec[2] };
  t8_vec_ax (vec_check, x);
  if (t8_vec_dist (vec_check, b) > tolerance) {
    /* Point does not lie on the line. */
    return 0;
  }
  /* The point is on the line. */
  return 1;
}

/**
 * Check if a point is inside of a triangle described by a point \a p_0 and two vectors \a v and \a w. 
 * 
 * \param[in] p_0         The first vertex of a triangle
 * \param[in] v           The vector from p_0 to p_1 (second vertex in the triangle)
 * \param[in] w           The vector from p_0 to p_2 (third vertex in the triangle)
 * \param[in] point       The coordinates of the point to check
 * \param[in] tolerance   A double > 0 defining the tolerance
 * \return                0 if the point is outside, 1 otherwise.  
 */
static int
t8_triangle_point_inside (const double p_0[3], const double v[3], const double w[3], const double point[3],
                          const double tolerance)
{
  /* A point p is inside the triangle that is spanned
   * by the point p_0 and vectors v and w if and only if the linear system
   * vx + wy = point - p_0
   * has a solution with 0 <= x,y and x + y <= 1.
   *
   * We check whether such a solution exists by computing
   * certain determinants of 2x2 submatrizes of the 3x3 matrix
   *
   *  | v w e_3 | with v = p_1 - p_0, w = p_2 - p_0, and e_3 = (0 0 1)^t (third unit vector)
   */

  T8_ASSERT (tolerance > 0); /* negative values and zero are not allowed */
  double b[3];
  /* b = point - p_0 */
  t8_vec_axpyz (p_0, point, b, -1);

  /* Let d = det (v w e_3) */
  const double det_vwe3 = v[0] * w[1] - v[1] * w[0];

  /* The system has a solution, we need to compute it and
   * check whether 0 <= x,y and x + y <= 1 */
  /* x = det (b w e_3) / d
   * y = det (v b e_3) / d
   */
  const double x = (b[0] * w[1] - b[1] * w[0]) / det_vwe3;
  const double y = (v[0] * b[1] - v[1] * b[0]) / det_vwe3;

  if (x < -tolerance || y < -tolerance || x + y > 1 + tolerance) {
    /* The solution is not admissible.
     * x < 0 or y < 0 or x + y > 1 */
    return 0;
  }
  /* The solution may be admissible, but we have to
   * check whether the result of
   *  (p_1 - p_0)x + (p_2 - p_0)y ( = vx + wy)
   * is actually p - p_0.
   * Since the system of equations is overrepresented (3 equations, 2 variables)
   * this may actually break.
   * If it breaks, it will break in the z coordinate of the result.
   */
  const double z = v[2] * x + w[2] * y;
  /* Must match the last coordinate of b = p - p_0 */
  if (fabs (z - b[2]) > tolerance) {
    /* Does not match. Point lies outside. */
    return 0;
  }
  /* All checks passed. Point lies inside. */
  return 1;
}

/** Check if a point lays on the inner side of a plane of a bilinearly interpolated volume element. 
 * the plane is described by a point and the normal of the face. 
 * \param[in] point_on_face   A point on the plane
 * \param[in] face_normal     The normal of the face
 * \param[in] point           The point to check
 * \return                    0 if the point is outside, 1 otherwise.                   
 */
static int
t8_plane_point_inside (const double point_on_face[3], const double face_normal[3], const double point[3])
{
  /* Set x = x - p */
  double pof[3] = { point_on_face[0], point_on_face[1], point_on_face[2] };
  t8_vec_axpy (point, pof, -1);
  /* Compute <x-p,n> */
  const double dot_product = t8_vec_dot (pof, face_normal);
  if (dot_product < 0) {
    /* The point is on the wrong side of the plane */
    return 0;
  }
  return 1;
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
                                                        const int num_points, int *is_inside, const double tolerance)
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
    for (int iface = 0; iface < num_faces; ++iface) {
      double face_normal[3];
      /* Compute the outer normal n of the face */
      t8_forest_element_face_normal (forest, ltreeid, element, iface, face_normal);
      /* Compute a point x on the face */
      const int afacecorner = ts->t8_element_get_face_corner (element, iface, 0);
      double point_on_face[3];
      t8_forest_element_coordinate (forest, ltreeid, element, afacecorner, point_on_face);
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        const int is_inside_iface = t8_plane_point_inside (point_on_face, face_normal, &points[ipoint * 3]);
        if (is_inside_iface == 0) {
          /* Point is on the outside of face iface. Update is_inside */
          is_inside[ipoint] = 0;
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
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_linear_new (int dimension)
{
  t8_geometry_linear *geom = new t8_geometry_linear (dimension);
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
