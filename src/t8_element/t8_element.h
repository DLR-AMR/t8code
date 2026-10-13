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
#include <t8_eclass/t8_eclass.h>
#include <t8_element/t8_element_shape.h>
#ifdef __cplusplus
#include <t8_types/t8_vec.hxx>
#endif

/** We want to export the whole implementation to be callable from "C". */
T8_EXTERN_C_BEGIN ();

/** Opaque structure for a generic element, only used as pointer.
 * Implementations are free to cast it to their internal data structure.
 */
typedef struct t8_element t8_element_t;

/* clang-format off */
/** Reference coordinates of each vertex of each element class. */
#define T8_ELEMENT_CORNER_REF_COORDS_VALUES {                                                     \
  { { 0, 0, 0 } },                                                                  /* T8_ECLASS_VERTEX */   \
  { { 0, 0, 0 }, { 1, 0, 0 } },                                                     /* T8_ECLASS_LINE */     \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 } },                           /* T8_ECLASS_QUAD */     \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } },                                        /* T8_ECLASS_TRIANGLE */ \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },                                                      \
    { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 } },                           /* T8_ECLASS_HEX */      \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, 1 } },                           /* T8_ECLASS_TET */      \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 } }, /* T8_ECLASS_PRISM */    \
  { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 1, 1 } }               /* T8_ECLASS_PYRAMID */  \
}

/** Reference coordinates of the centroid of each element class. */
#define T8_ELEMENT_CENTROID_REF_COORDS_VALUES {      \
  { 0, 0, 0 },               /* T8_ECLASS_VERTEX */   \
  { 0.5, 0, 0 },             /* T8_ECLASS_LINE */     \
  { 0.5, 0.5, 0 },           /* T8_ECLASS_QUAD */     \
  { 2. / 3., 1. / 3., 0 },   /* T8_ECLASS_TRIANGLE */ \
  { 0.5, 0.5, 0.5 },         /* T8_ECLASS_HEX */      \
  { 0.75, 0.25, 0.5 },       /* T8_ECLASS_TET */      \
  { 2. / 3., 1. / 3., 0.5 }, /* T8_ECLASS_PRISM */    \
  { 0.6, 0.6, 0.2 }          /* T8_ECLASS_PYRAMID */  \
}
/* clang-format on */

#ifdef __cplusplus
/* constexpr variables for cpp. As in t8_eclass.h, they are wrapped in a namespace to have a different symbol
 * than the C variables. T8_EXTERN_C is disabled, because it disables the namespace. */

T8_EXTERN_C_END ();

namespace t8cpp
{
/** This array holds the reference coordinates of each vertex of each element.
 *  It can e.g. be used with the \ref t8_scheme::element_get_reference_coords function.
 *  Usage: t8_element_corner_ref_coords[eclass][vertex]
 */
inline constexpr t8_3D_vec t8_element_corner_ref_coords[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS]
  = T8_ELEMENT_CORNER_REF_COORDS_VALUES;

/** This array holds the reference coordinates of the centroid of each element.
 *  It can e.g. be used with the \ref t8_scheme::element_get_reference_coords function.
 *  Usage: t8_element_centroid_ref_coords[eclass]
 */
inline constexpr t8_3D_vec t8_element_centroid_ref_coords[T8_ECLASS_COUNT] = T8_ELEMENT_CENTROID_REF_COORDS_VALUES;
} /* namespace t8cpp */

using namespace t8cpp;

T8_EXTERN_C_BEGIN ();

#else /* !__cplusplus */

/** This array holds the reference coordinates of each vertex of each element.
 *  It can e.g. be used with the \ref t8_element_get_reference_coords function.
 *  Usage: t8_element_corner_ref_coords[eclass][vertex][dimension]
 */
extern const double t8_element_corner_ref_coords[T8_ECLASS_COUNT][T8_ECLASS_MAX_CORNERS][3];

/** This array holds the reference coordinates of the centroid of each element.
 *  It can e.g. be used with the \ref t8_element_get_reference_coords function.
 *  Usage: t8_element_centroid_ref_coords[eclass][dimension]
 */
extern const double t8_element_centroid_ref_coords[T8_ECLASS_COUNT][3];

#endif /* !__cplusplus */

/* Undefine values so that they do not leak into other files.
 * They can be kept using KEEP_ELEMENT_VALUE_DEFINITIONS. This is needed
 * to use them in the t8_element.c file. */
#ifndef KEEP_ELEMENT_VALUE_DEFINITIONS
#undef T8_ELEMENT_CORNER_REF_COORDS_VALUES
#undef T8_ELEMENT_CENTROID_REF_COORDS_VALUES
#endif

/** End of code that is callable from "C".*/
T8_EXTERN_C_END ();

#endif /* !T8_ELEMENT_H */
