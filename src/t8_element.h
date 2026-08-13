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

/** \file t8_element.h
 * This file defines the opaque element structure and provides some
 * constants for element classes.
 */

#ifndef T8_ELEMENT_H
#define T8_ELEMENT_H

#include <t8.h>
#include <t8_eclass.h>
#include <t8_element_shape.h>

T8_EXTERN_C_BEGIN ();

/** Opaque structure for a generic element, only used as pointer.
 * Implementations are free to cast it to their internal data structure.
 */
typedef struct t8_element t8_element_t;

/** This array holds the reference coordinates of each vertex of each element.
 *  It can e.g. be used with the \ref t8_scheme::element_get_reference_coords function.
 *  Usage: t8_element_corner_ref_coords[eclass][vertex][dimension]
 */
extern const double t8_element_corner_ref_coords[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS][3];

/** This array holds the reference coordinates of the centroid of each element.
 *  It can e.g. be used with the \ref t8_scheme::element_get_reference_coords function.
 *  Usage: t8_element_centroid_ref_coords[eclass][dimension]
 */
extern const double t8_element_centroid_ref_coords[T8_ECLASS_COUNT][3];

T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
