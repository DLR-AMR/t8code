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

#ifndef T8_DVERTEX_H
#define T8_DVERTEX_H

/** \file t8_dvertex.h
 * TODO: document this.
 */

#include <t8.h>

/** The number of children that a vertex is refined into. */
#define T8_DVERTEX_CHILDREN 1

/** The number of faces of a vertex. */
#define T8_DVERTEX_FACES 0

/** The number of face children of a vertex. */
#define T8_DVERTEX_FACE_CHILDREN 0

/** The length of a vertex root tree */
#define T8_DVERTEX_ROOT_LEN 0

/** The maximum refinement level allowed for a vertex.
 * The max level is lower than 255 so that we can use an uint8_t
 * to iterate to maxlevel:
 * for (uint8_t level = 0; level <= T8_DVERTEX_MAXLEVEL; ++level)
 * Otherwise, uint8_t would overflow after 255 and we would have an infinite loop.
*/
#define T8_DVERTEX_MAXLEVEL 254

/** The data container describing a refined element in a refined tree for the vertex element class. */
typedef struct t8_dvertex
{
  uint8_t level;
} t8_dvertex_t;

#endif /* T8_DVERTEX_H */
