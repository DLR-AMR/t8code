/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/** \file t8_type.h
 */

#ifndef T8_TYPE_H
#define T8_TYPE_H

#include <t8.h>

/** This enumeration contains all possible element types. */
typedef enum t8_type
{
  T8_TYPE_FIRST = 0,
  T8_TYPE_VERTEX = T8_TYPE_FIRST,
  T8_TYPE_LINE,
  T8_TYPE_QUAD,
  T8_TYPE_TRIANGLE,
  T8_TYPE_HEX,
  T8_TYPE_TET,
  T8_TYPE_PRISM,
  T8_TYPE_PYRAMID,
  T8_TYPE_LAST
}
t8_type_t;

/** Map each of the element types to its dimension. */
extern const int    t8_type_to_dimension[T8_TYPE_LAST];

/** For each of the element types, count the boundary points by type. */
extern const int    t8_type_boundary_count[T8_TYPE_LAST][T8_TYPE_LAST];

/** Query the type and count of boundary points.
 * \param [in] thetype          We query a point of this type.
 * \param [in] min_dim          Ignore boundary points of lesser dimension.
 *                              The ignered points get a count value of 0.
 * \param [out] per_type        Array of length T8_TYPE_LAST to be filled
 *                              with the count of the boundary objects,
 *                              counted per each of the element types.
 * \return                      The count over all boundary points.
 */
int                 t8_type_count_boundary (t8_type_t thetype,
                                            int min_dim, int *per_type);

#endif /* !T8_ELEMENT_H */
