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

/** \file t8_dhex.h */

#ifndef T8_DHEX_H
#define T8_DHEX_H

#include <t8.h>

/** Spatial dimension of a hexahedron. */
#define T8_DHEX_DIM 3

/** The number of children that a hexahedron is refined into. */
#define T8_DHEX_CHILDREN 8

/** The maximum refinement level allowed for a node within a hexahedron. */
#define T8_DHEX_MAXLEVEL 19

/** The maximum refinement level allowed for a hexahedron. */
#define T8_DHEX_QMAXLEVEL 18

/** The number of faces of a hexahedron. */
#define T8_DHEX_FACES (2 * T8_DHEX_DIM)

/** The length of the root hexahedron in integer coordinates. */
#define T8_DHEX_ROOT_LEN ((t8_qcoord_t) 1 << T8_DHEX_MAXLEVEL)

/** The length of a hexahedron at a given level in integer coordinates. */
#define T8_DHEX_LEN(l) ((t8_qcoord_t) 1 << (T8_DHEX_MAXLEVEL - (l)))

/** The offset of the highest (farthest from the origin) hexahedron at level l. */
#define T8_DHEX_LAST_OFFSET(l) (T8_DHEX_ROOT_LEN - T8_DHEX_LEN (l))

/** This data type encodes a hexahedron. */
typedef struct t8_dhex
{
  t8_qcoord_t         x, y, z;
  int8_t              level;
}
t8_dhex_t;

#endif /* T8_DHEX_H */
