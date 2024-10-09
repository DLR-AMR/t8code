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

#include <t8_schemes/t8_default/t8_default_hex/t8_dhex_bits.h>
#include <p8est_bits.h>

void
t8_dhex_compute_reference_coords (const t8_dhex_t *elem, const double *ref_coords, const size_t num_coords,
                                  const size_t padding, double *out_coords)
{
  const p8est_quadrant_t *q1 = (const p8est_quadrant_t *) elem;

  /* Get the length of the quadrant */
  const p4est_qcoord_t len = P8EST_QUADRANT_LEN (q1->level);

  for (size_t coord = 0; coord < num_coords; ++coord) {
    const size_t offset_ref = (t8_eclass_to_dimension[T8_ECLASS_HEX] + padding) * coord;
    const size_t offset_out = t8_eclass_to_dimension[T8_ECLASS_HEX] * coord;
    /* Compute the x, y and z coordinates of the point depending on the
     * reference coordinates */
    out_coords[offset_out + 0] = q1->x + ref_coords[offset_ref + 0] * len;
    out_coords[offset_out + 1] = q1->y + ref_coords[offset_ref + 1] * len;
    out_coords[offset_out + 2] = q1->z + ref_coords[offset_ref + 2] * len;

    /* We divide the integer coordinates by the root length of the hex
     * to obtain the reference coordinates. */
    out_coords[offset_out + 0] /= (double) P8EST_ROOT_LEN;
    out_coords[offset_out + 1] /= (double) P8EST_ROOT_LEN;
    out_coords[offset_out + 2] /= (double) P8EST_ROOT_LEN;
  }
}
