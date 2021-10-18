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

void                t8_geom_compute_linear_geometry (t8_eclass_t tree_class,
                                                     const double
                                                     *tree_vertices,
                                                     const double *ref_coords,
                                                     double out_coords[3]);

/** Given function values at the four edge points of a unit square and
 * a point within that square, interpolate the function value at this point.
 * \param [in]    vertex  An array of size at least dim giving the coordinates of the vertex to interpolate
 * \param [in]    corner_values An array of size 2^dim * 3, giving for each corner (in zorder) of
 *                        the unit square/cube its function values in 3D space.
 * \param [out]   evaluated_function An array of size 3, on output the function values
 *                        at \a vertex are stored here.
 */
void                t8_geom_bilinear_interpolation (const double *vertex,
                                                    const double
                                                    *corner_values, int dim,
                                                    double
                                                    *evaluated_function);

T8_EXTERN_C_END ();

#endif /* !T8_GEOMETRY_HELPERS_H! */
