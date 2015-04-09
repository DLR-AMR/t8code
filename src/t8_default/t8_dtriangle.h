/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#ifndef T8_DTRIANGLE_H
#define T8_DTRIANGLE_H

/** The number of children that a triangle is refined into. */
#define T8_DTRIANGLE_CHILDREN 8

/** The number of faces of a triangle. */
#define T8_DTRIANGLE_FACES 4

/** The maximum refinement level allowed for a triangle. */
#define T8_DTRIANGLE_MAXLEVEL 30

/** The length of the root triangle in integer coordinates. */
#define T8_DTRIANGLE_ROOT_LEN (1 << (T8_DTET_MAXLEVEL))

/** The length of a triangle at a given level in integer coordinates. */
#define T8_DTRIANGLE_LEN(l) (1 << (T8_DTRIANGLE_MAXLEVEL - (l)))
typedef int32_t     t8_dtriangle_coord_t;

typedef struct t8_dtriangle
{
  int8_t              level;
  /* add coordinates etc. here */
  int8_t              type;
  t8_dtriangle_coord_t x, y;
}
t8_dtriangle_t;

typedef int8_t      t8_dtriangle_type_t;

#endif /* T8_DTRIANGLE_H */
