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

#ifndef T8_DHEX_H
#define T8_DHEX_H

/** \file t8_dhex.h
 * TODO: document this.
 */

#include <t8.h>

/** The number of children that a hex is refined into. */
#define T8_DHEX_CHILDREN 8

/** The number of faces of a hex. */
#define T8_DHEX_FACES 6

/** The number of vertices of a hex. */
#define T8_DHEX_VERTICES 8

/** The number of children at a face of a hex. */
#define T8_DHEX_FACE_CHILDREN 4

/** The maximum refinement level allowed for a hex. */
#define T8_DHEX_MAXLEVEL 21

/** The length of the root hex in integer coordinates. */
#define T8_DHEX_ROOT_LEN (1 << (T8_DHEX_MAXLEVEL))

/** The length of a hex at a given level in integer coordinates. */
#define T8_DHEX_LEN(l) (1 << (T8_DHEX_MAXLEVEL - (l)))

/** Type of an integer coordinate for a node of a hex element. */
typedef int32_t t8_dhex_coord_t;

/** The data container describing a refined element in a refined tree for the hex element class. */
typedef struct t8_dhex
{
  int8_t level;
  t8_dhex_coord_t x; /**< The x integer coordinate of the anchor node. */
  t8_dhex_coord_t y; /**< The y integer coordinate of the anchor node. */
  t8_dhex_coord_t z; /**< The z integer coordinate of the anchor node. */
} t8_dhex_t;

#endif /* T8_DHEX_H */
