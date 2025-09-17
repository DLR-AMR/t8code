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

#ifndef T8_DTRI_H
#define T8_DTRI_H

/** \file t8_dtri.h
 * TODO: document this.
 */

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** The number of children that a triangle is refined into. */
#define T8_DTRI_CHILDREN 4

/** The number of faces of a triangle. */
#define T8_DTRI_FACES 3

/** The number of children that a face of a triangle is refined to. */
#define T8_DTRI_FACE_CHILDREN 2

/** The number of corners of a triangle */
#define T8_DTRI_CORNERS 3

/** The maximum refinement level allowed for a triangle. */
#define T8_DTRI_MAXLEVEL 29

/** The length of the root triangle in integer coordinates. */
#define T8_DTRI_ROOT_LEN (1 << (T8_DTRI_MAXLEVEL))

/** The length of a triangle at a given level in integer coordinates. */
#define T8_DTRI_LEN(l) (1 << (T8_DTRI_MAXLEVEL - (l)))

/** The number of different types of triangles. */
#define T8_DTRI_NUM_TYPES 2

/** The length of a line divided by the length of a triangle.
 *  This is useful to convert boundary coordinates from tri to line. */
#define T8_DLINE_ROOT_BY_DTRI_ROOT (1 << (T8_DLINE_MAXLEVEL - T8_DTRI_MAXLEVEL))

/** Type for the (integer) type of a triangular element. */
typedef int8_t t8_dtri_type_t;
/** Type of an integer coordinate for a node of a triangular element. */
typedef int32_t t8_dtri_coord_t;

/** The data container describing a refined element in a refined tree for the triangular element class. */
typedef struct t8_dtri
{
  int8_t level;        /**< The refinement level of the element relative to the root at level 0. */
  t8_dtri_type_t type; /**< Type of the triangle (0 or 1). */
  t8_dtri_coord_t x;   /**< The x integer coordinate of the anchor node. */
  t8_dtri_coord_t y;   /**< The y integer coordinate of the anchor node. */
} t8_dtri_t;

T8_EXTERN_C_END ();

#endif /* T8_DTRI_H */
