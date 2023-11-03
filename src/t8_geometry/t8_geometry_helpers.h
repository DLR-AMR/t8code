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

/** \file t8_geometry_helpers.h
 * Defines t8code internal functions that are useful for geometry
 * implementations.
 */

#ifndef T8_GEOMETRY_HELPERS_H
#define T8_GEOMETRY_HELPERS_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Compute the linear geometry of a tree at a given reference coordinate.
 * \param [in]    tree_class     The eclass of the tree.
 * \param [in]    tree_vertices  Array with the tree vertex coordinates.
 * \param [in]    ref_coords     The reference coordinates of the point.
 * \param [in]    num_coords     Number of points to evaluate.
 * \param [out]   out_coords     The output coordinates.
 */
void
t8_geom_compute_linear_geometry (t8_eclass_t tree_class, const double *tree_vertices, const double *ref_coords,
                                 const size_t num_coords, double *out_coords);

/** Compute the linear, axis-aligned geometry of a tree at a given reference coordinate.
 *  This function is faster than \ref t8_geom_compute_linear_geometry, but only works
 *  for axis-aligned trees of \ref T8_ECLASS_LINE, \ref T8_ECLASS_QUAD and \ref T8_ECLASS_HEX.
 * \param [in]    tree_class     The eclass of the tree.
 * \param [in]    tree_vertices  Array with the tree vertex coordinates.
 * \param [in]    ref_coords     The reference coordinates of the point.
 * \param [in]    num_coords     Number of points to evaluate.
 * \param [out]   out_coords     The output coordinates.
 */
void
t8_geom_compute_linear_axis_aligned_geometry (t8_eclass_t tree_class, const double *tree_vertices,
                                              const double *ref_coords, const size_t num_coords, double *out_coords);

/** Interpolates linearly between 2, bilinearly between 4 or trilineraly between 8 points.
 * \param [in]    coefficients        An array of size at least dim giving the coefficients used for the interpolation
 * \param [in]    corner_values       An array of size 2^dim * 3, giving for each corner (in zorder) of
 *                                    the unit square/cube its function values in space.
 * \param [in]    corner_value_dim    The dimension of the \a corner_values.
 * \param [in]    interpolation_dim   The dimension of the interpolation (1 for linear, 2 for bilinear, 3 for trilinear)
 * \param [out]   evaluated_function  An array of size \a corner_value_dim, on output the result of the interpolation.
 */
void
t8_geom_linear_interpolation (const double *coefficients, const double *corner_values, int corner_value_dim,
                              int interpolation_dim, double *evaluated_function);

/** Triangular interpolation between 3 points (triangle) or 4 points (tetrahedron) using barycentric coordinates.
 * \param [in]    coefficients        An array of size \a interpolation_dim giving the coefficients used for the interpolation
 * \param [in]    corner_values       An array of size 
                                                       3 * \a corner_value_dim for \a interpolation_dim == 2 or
                                                       4 * \a corner_value_dim for \a interpolation_dim == 3, 
                                                       giving the function values of the triangle/tetrahedron for each corner (in zorder)
 * \param [in]    corner_value_dim    The dimension of the \a corner_values.
 * \param [in]    interpolation_dim   The dimension of the interpolation (2 for triangle, 3 for tetrahedron)
 * \param [out]   evaluated_function  An array of size \a corner_value_dim, on output the result of the interpolation.
 */
void
t8_geom_triangular_interpolation (const double *coefficients, const double *corner_values, int corner_value_dim,
                                  int interpolation_dim, double *evaluated_function);

/** Copies the vertex coordinates of a tree face in zorder into a separate array.
 * \param [in]    tree_class     The eclass of the tree.
 * \param [in]    tree_vertices  Array with the tree vertex coordinates.
 * \param [in]    face_index     Index of the face, which vertices should be copied.
 * \param [in]    dim            The dimension of the face vertices.
 * \param [out]   face_vertices  Coordinates of the face vertices in zorder.
 */
void
t8_geom_get_face_vertices (t8_eclass_t tree_class, const double *tree_vertices, int face_index, int dim,
                           double *face_vertices);

/** Copies the vertex coordinates of a tree edge in zorder into a separate array.
 * \param [in]    tree_class     The eclass of the tree.
 * \param [in]    tree_vertices  Array with the tree vertex coordinates.
 * \param [in]    edge_index     Index of the edge, which vertices should be copied.
 * \param [in]    dim            The dimension of the edge vertices.
 * \param [out]   edge_vertices  Coordinates of the edge vertices in zorder.
 */
void
t8_geom_get_edge_vertices (t8_eclass_t tree_class, const double *tree_vertices, int edge_index, int dim,
                           double *edge_vertices);

/** Calculates a point of intersection in a triangular reference space. The intersection is the extension
 * of a straight line passing through a reference point and the opposite vertex of the edge.
 *     /|\
 *    / | \     o -> reference point
 *   /  o  \    x -> intersection point
 *  /   |   \
 * /____x____\
 * \param [in]    edge_index        Index of the edge, the intersection lies on.
 * \param [in]    ref_coords        Array containing the coordinates of the reference point.
 * \param [out]   ref_intersection  Coordinates of the intersection point.
 */
void
t8_geom_get_ref_intersection (int edge_index, const double *ref_coords, double ref_intersection[2]);

/** Calculates the scaling factor for edge displacement along a triangular tree face
 *  depending on the position of the global reference point.
 * \param [in]    edge_index          Index of the edge, whose displacement should be scaled.
 * \param [in]    tree_vertices       Array with the tree vertex coordinates.
 * \param [in]    glob_intersection   Array containing the coordinates of the intersection point
 *                                    of a line drawn from the opposite vertex through the
 *                                    glob_ref_point onto the edge with edge_index.
 * \param [in]    glob_ref_point      Array containing the coordinates of the reference point
 *                                    mapped into the global space.
 */
double
t8_geom_get_triangle_scaling_factor (int edge_index, const double *tree_vertices, const double *glob_intersection,
                                     const double *glob_ref_point);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_HELPERS_H! */
