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

#ifndef T8_TRANSFORMATION_MATRIX
#define T8_TRANSFORMATION_MATRIX

#include <sc_dmatrix.h>

/**
 * We define several projection matrices to transform coordinates from the world
 * coordinate system into a camera coordinate-system and to project the onto a plane
 * As the plane hase x-y-axis it is convenient if the camera-system coincides with
 * this coordinate system. Therfore the z-axis is facing "backwards" and the y-axis
 * is faceing forward. 
 * 
 * From a cameras-point of view:
 *         y               /| 
 *         ^              / |
 *         |             /  |           
 *         |            |   |
 * z<------|           y|  /
 *        /             | /x
 *       /              |/
 *      /          projection-plane
 *      x 
 */

/**
 * We use a similar (slightly simplified) camera-modell, where the up-vector is 
 * assumend to always be the z-axis of the world-coordinate-system (the forest-coordinate system)
 * 
 * \param[in] cam       The coordinates of the camera in the world coordinate system
 * \param[in] ref_point The point the camera is focussing
 * \param[in] cam_y     The up direction in world-coordinates 
 * \param[in, out] transformation A zero-initialized matrix on input. On output the inverse of the camera-transformation matrix, mapping points from
 *                      the world-coordinate-system into the camera-coordinate system.
 */
void                inverse_camera_transformation (const double cam[3],
                                                   const double ref_point[3],
                                                   const double up[3],
                                                   double
                                                   transformation[4][4]);

/**
 * Compute the perspective projection matrix, mapping a a view-volume into the cube [0,1]^3. 
 * 
 * \param[in] width width of the image plane
 * \param[in] height height of the image plane
 * \param[in] near distance of the image plane along the z-axis
 * \param[in] far distance of the furthest displayed area along the z-axis
 * \param[in, out] projection A zero-initialized matrix. On output the perspective projection matrix, projecting points in the
 *                  area defined by the image-plane, near and far-field into the cube [-1, 1]^3.
 *                  a pinhole camera model. A parallel projection can than be used to project the points
 *                  down to a plane. 
 * 
 *                                                 /|
 *                                                / |
 *         y               /|                    /  |
 *         ^              / |height             /   |
 *         |             /  |                   |   |
 *         |            |   |                   |   |
 * z<------|------near--|  /                    |   /
 *        /             | /width                |  /
 *       /--------------|/----------far---------| /
 *      /                                       |/
 *      x 
 * 
 * 
 */
void                perspective_projection (const double width,
                                            const double height,
                                            const double near,
                                            const double far,
                                            double projection[4][4]);

/**
 * A matrix-vector product for projection matrices and 3D vectors. 
 *
 * \param[in] mat A projection matrix
 * \param[in] vec A 3D vector
 * \param[in, out] out_vec On output the matrix-vector product of \a mat * \a vec. 
 */
void                mat4d_vec_multi (const double mat[4][4],
                                     const double vec[3], double out_vec[3]);
#endif /* T8_TRANSFORMATION_MATRIX */
