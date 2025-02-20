/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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

/** \file t8_vec.h
 * Provide vector operations for 3D vectors.
 */

#ifndef T8_VEC_H
#define T8_VEC_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Vector norm.
 * \param [in] vec  A 3D vector.
 * \return          The norm of \a vec.
 */
double
t8_norm (const double vec[3]);

/** Normalize a vector.
 * \param [in,out] vec  A 3D vector.
 */
void
t8_normalize (double vec[3]);

/** Make a copy of a vector.
 * \param [in]  vec_in
 * \param [out] vec_out
 */
void
t8_copy (const double vec_in[3], double vec_out[3]);

/** Euclidean distance of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The euclidean distance.
 *                     Equivalent to norm (X-Y).
 */
double
t8_dist (const double vec_x[3], const double vec_y[3]);

/** Compute X = alpha * X
 * \param [in,out] vec_x  A 3D vector. On output set to \a alpha * \a vec_x.
 * \param [in]     alpha  A factor.
 */
void
t8_ax_c_interface (double vec_x[3], const double alpha);

/** Compute Y = alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_z  On output set to \a alpha * \a vec_x.
 * \param [in]  alpha  A factor.
 */
void
t8_axy_c_interface (const double vec_x[3], double vec_y[3], const double alpha);

/** Y = alpha * X + b
 * \param [in]  vec_x  A 3D vector.
 * \param [out] vec_y  On input, a 3D vector.
 *                     On output set to \a alpha * \a vec_x + \a b.
 * \param [in]  alpha  A factor.
 * \param [in]  b      An offset.
 * \note It is possible that vec_x = vec_y on input to overwrite x
 */
void
t8_axb_c_interface (const double vec_x[3], double vec_y[3], const double alpha, const double b);
/** Y = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in,out] vec_y On input, a 3D vector.
 *                      On output set \a to vec_y + \a alpha * \a vec_x
 * \param [in]  alpha  A factor.
 */
void
t8_axpy_c_interface (const double vec_x[3], double vec_y[3], const double alpha);

/** Z = Y + alpha * X
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] vec_z  On output set \a to vec_y + \a alpha * \a vec_x
 */
void
t8_axpyz_c_interface (const double vec_x[3], const double vec_y[3], double vec_z[3], const double alpha);

/** Dot product of X and Y.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \return             The dot product \a vec_x * \a vec_y
 */
double
t8_dot_c_interface (const double vec_x[3], const double vec_y[3]);

/** Cross product of X and Y
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
void
t8_cross_3D_c_interface (const double vec_x[3], const double vec_y[3], double cross[3]);

/** Cross product of X and Y
 * \param [in]  vec_x  A 2D vector.
 * \param [in]  vec_y  A 2D vector.
 * \param [out] cross  On output, the cross product of \a vec_x and \a vec_y.
 */
double
t8_cross_2D_c_interface (const double vec_x[2], const double vec_y[2]);

/** Compute the difference of two vectors.
 * \param [in]  vec_x  A 3D vector.
 * \param [in]  vec_y  A 3D vector.
 * \param [out] diff   On output, the difference of \a vec_x and \a vec_y.
 */
void
t8_diff_c_interface (const double vec_x[3], const double vec_y[3], double diff[3]);

/**
 * Check the equality of two vectors elementwise 
 * 
 * \param[in] vec_x 
 * \param[in] vec_y 
 * \param[in] tol 
 * \return true, if the vectors are equal up to \a tol 
 */
int
t8_eq_c_interface (const double vec_x[3], const double vec_y[3], const double tol);

/** Rescale a vector to a new length.
 * \param [in,out] vec  A 3D vector.
 * \param [in]  new_length  New length of the vector.
 */
void
t8_rescale_c_interface (double vec[3], const double new_length);

/** Compute the normal of a triangle given by its three vertices.
 * \param [in]  p1  A 3D vector.
 * \param [in]  p2  A 3D vector.
 * \param [in]  p3  A 3D vector.
 * \param [out] Normal vector of the triangle. (Not necessarily of length 1!)
 */
void
t8_normal_of_tri_c_interface (const double p1[3], const double p2[3], const double p3[3], double normal[3]);

/** Compute an orthogonal coordinate system from a given vector.
 * \param [in]   v1 3D vector.
 * \param [out]  v2 3D vector.
 * \param [out]  v3 3D vector.
 */
void
t8_orthogonal_tripod_c_interface (const double v1[3], double v2[3], double v3[3]);
/** Swap the components of two vectors.
 * \param [in,out]  p1  A 3D vector.
 * \param [in,out]  p2  A 3D vector.
 */
void
t8_swap_c_interface (double p1[3], double p2[3]);

T8_EXTERN_C_END ();

#endif /* T8_VEC_H */
