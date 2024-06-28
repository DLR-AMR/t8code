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

/** \file t8_dquad_bits.h
 * Definitions of quad-specific functions.
 */

#ifndef T8_DQUAD_BITS_H
#define T8_DQUAD_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad.h>

T8_EXTERN_C_BEGIN ();

/** Convert points in the reference space of a quad element to points in the
 *  reference space of the tree (level 0) embedded in \f$ [0,1]^2 \f$.
 * \param [in]  elem       Input quad.
 * \param [in]  ref_coords The reference coordinates in the quad
 *                         (\a num_coords times \f$ [0,1]^2 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [out] out_coords An array of \a num_coords x 2 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the quad.
 */
void
t8_dquad_compute_reference_coords (const p4est_quadrant_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords);

/** Compute the childid-th child in Morton order of a quad.
 * \param [in] elem       Input quad.
 * \param [in] childid    The id of the child, in 0 - 3, in Morton order.
 * \param [in,out] child  Existing quad whose data will be filled with the data of elem's childid-th child.
 */
void
t8_dquad_child (const p4est_quadrant_t *elem, int childid, p4est_quadrant_t *child);

/** Compute the parent of a quad.
 * \param [in]  elem Input quad.
 * \param [in,out] parent Existing quad whose data will
 *                  be filled with the data of elem's parent.
 * \note \a elem may point to the same quad as \a parent.
 */
void
t8_dquad_parent (const p4est_quadrant_t *elem, p4est_quadrant_t *parent);

/** Compute the sibid-th sibling in Morton order of a quad.
 * \param [in] elem         Input quad.
 * \param [in] sibid        The id of the sibling, in 0 - 3, in Morton order.
 * \param [in,out] sibling  Existing quad whose data will be filled with the data of elem's sibid-th sibling.
 */
void
t8_dquad_sibling (const p4est_quadrant_t *elem, int sibid, p4est_quadrant_t *sibling);

/** Copy all values from one quad to another.
 * \param [in] elem     The quad to be copied.
 * \param [in,out] dest Existing quad whose data will be filled with the data
 *                   of \a elem.
 */
void
t8_dquad_copy (const p4est_quadrant_t *source, p4est_quadrant_t *dest);

/** Computes the successor of a quad in a uniform grid of level \a level.
 * \param [in] elem      Quad whose id will be computed.
 * \param [in,out] succ  Existing quad whose data will be filled with the
 *                       data of \a elem's successor on level \a level.
 * \param [in] level     Level of uniform grid to be considered.
 * \param [in] multilevel If \a multilevel is 1, ancestors can also be successors.
 */
void
t8_dquad_successor (const p4est_quadrant_t *elem, p4est_quadrant_t *succ, const int level, const int multilevel);

T8_EXTERN_C_END ();

#endif /* T8_DQUAD_BITS_H */
