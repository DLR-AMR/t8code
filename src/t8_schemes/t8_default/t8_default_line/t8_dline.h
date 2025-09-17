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

#ifndef T8_DLINE_H
#define T8_DLINE_H

/** \file t8_dline.h
 * TODO: document this.
 */

#include <t8.h>

/** The number of children that a line is refined into. */
#define T8_DLINE_CHILDREN 2

/** The number of faces of a line. */
#define T8_DLINE_FACES 2

/** The number of children at a face of a line. */
#define T8_DLINE_FACE_CHILDREN 1

/** The maximum refinement level allowed for a line. */
#define T8_DLINE_MAXLEVEL 30

/** The length of the root line in integer coordinates. */
#define T8_DLINE_ROOT_LEN (1 << (T8_DLINE_MAXLEVEL))

/** The length of a line at a given level in integer coordinates. */
#define T8_DLINE_LEN(l) (1 << (T8_DLINE_MAXLEVEL - (l)))

/** Type of an integer coordinate for a node of a line element. */
typedef int32_t t8_dline_coord_t;

/** The data container describing a refined element in a refined tree for the line element class. */
typedef struct t8_dline
{
  t8_dline_coord_t x; /**< The integer coordinate of the anchor node. */
  int8_t level;       /**< The refinement level of the element relative to the root at level 0. */
} t8_dline_t;

#endif /* T8_DLINE_H */
