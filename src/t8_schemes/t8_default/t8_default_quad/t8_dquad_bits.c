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

#include <p4est_bits.h>
#include <t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.h>
#include <p4est_bits.h>

void
t8_dquad_compute_reference_coords (const p4est_quadrant_t *elem, const double *ref_coords, const size_t num_coords,
                                   double *out_coords)
{
  const p4est_qcoord_t h = P4EST_QUADRANT_LEN (elem->level);

  for (size_t icoord = 0; icoord < num_coords; ++icoord) {
    const size_t offset_2d = icoord * 2;
    const size_t offset_3d = icoord * 3;
    out_coords[offset_2d + 0] = elem->x + ref_coords[offset_3d + 0] * h;
    out_coords[offset_2d + 1] = elem->y + ref_coords[offset_3d + 1] * h;

    out_coords[offset_2d + 0] /= (double) P4EST_ROOT_LEN;
    out_coords[offset_2d + 1] /= (double) P4EST_ROOT_LEN;
  }
}
