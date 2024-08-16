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

/** \file t8_dhex_bits.h
 * Definitions of hex-specific functions.
 */

#ifndef T8_DHEX_BITS_H
#define T8_DHEX_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_hex/t8_dhex.h>

T8_EXTERN_C_BEGIN ();

/** Convert points in the reference space of a hex element to points in the
 *  reference space of the tree (level 0) embedded in \f$ [0,1]^3 \f$.
 * \param [in]  elem       Input hex.
 * \param [in]  ref_coords The reference coordinates in the hex
 *                         (\a num_coords times \f$ [0,1]^3 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [in]  padding    Only used if \a num_coords > 1.
 *                         The amount of padding in \a ref_coords between array elements:
 *                         For elem dim 2 and input {x1, y1, x2, y2, ...} padding is 0.
 *                         For elem dim 2 and input {x1, y1, z1, x2, y2, z2, ...} padding is 1.
 *                         For elem dim 3 and input {x1, y1, z1, x2, y2, z2, ...} padding is 0.
 * \param [out] out_coords An array of \a num_coords x 3 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the hex.
 */
void
t8_dhex_compute_reference_coords (const t8_dhex_t *elem, const double *ref_coords, const size_t num_coords,
                                  const size_t padding, double *out_coords);

T8_EXTERN_C_END ();

#endif /* T8_DHEX_BITS_H */
