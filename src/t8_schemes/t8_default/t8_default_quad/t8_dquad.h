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

/** \file t8_dquad.h */

#ifndef T8_DQUAD_H
#define T8_DQUAD_H

#include <t8.h>

/** Spatial dimension of a quadrant. */
#define T8_DQUAD_DIM 2

/** The number of children that a quadrant is refined into. */
#define T8_DQUAD_CHILDREN 4

/** The maximum refinement level allowed for a node within a quadrant. */
#define T8_DQUAD_MAXLEVEL 30

/** The maximum refinement level allowed for a quadrant. */
#define T8_DQUAD_QMAXLEVEL 29

/** The number of faces of a quadrant. */
#define T8_DQUAD_FACES (2 * T8_DQUAD_DIM)

/** The length of the root quadrant in integer coordinates. */
#define T8_DQUAD_ROOT_LEN ((t8_qcoord_t) 1 << T8_DQUAD_MAXLEVEL)

/** The length of a quadrant at a given level in integer coordinates. */
#define T8_DQUAD_LEN(l) ((t8_qcoord_t) 1 << (T8_DQUAD_MAXLEVEL - (l)))

/** The offset of the highest (farthest from the origin) quadrant at level l. */
#define T8_DQUAD_LAST_OFFSET(l) (T8_DQUAD_ROOT_LEN - T8_DQUAD_LEN (l))

/** This data type encodes a quadrant. */
typedef struct t8_dquad
{
  t8_qcoord_t         x, y;
  int8_t              level;
}
t8_dquad_t;

#endif /* T8_DQUAD_H */
