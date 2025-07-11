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

/** \file t8_default_quad_bits.hxx
 * Definitions of quad-specific functions.
 */

#ifndef T8_DQUAD_BITS_H
#define T8_DQUAD_BITS_H

#include <t8_element.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad.hxx>

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
t8_dquad_compute_reference_coords (const t8_pquad_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords);

#endif /* T8_DQUAD_BITS_H */
