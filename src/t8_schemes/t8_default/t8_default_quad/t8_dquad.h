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

#ifndef T8_DQUAD_H
#define T8_DQUAD_H

/** \file t8_dhilbert.h
 * TODO: document this.
 */

#include <t8.h>

/** The number of children that a hilbert is refined into. */
#define T8_HILBERT_CHILDREN 4
#define T8_HILBERT_DIM 2

/** The maximum refinement level allowed for a hilbert. */
#define T8_HILBERT_MAXLEVEL 29

/** The length of the root hilbert in integer coordinates. */
#define T8_HILBERT_ROOT_LEN (1 << (T8_HILBERT_MAXLEVEL))

/** The length of a hilbert at a given level in integer coordinates. */
#define T8_HILBERT_LEN(l) (1 << (T8_HILBERT_MAXLEVEL - (l)))

typedef int32_t t8_hilbert_coord_t;
typedef int8_t t8_hilbert_type_t;

typedef struct t8_hilbert
{
  t8_hilbert_coord_t x; /**< The x integer coordinate of the anchor node. */
  t8_hilbert_coord_t y; /**< The y integer coordinate of the anchor node. */
  int8_t level;
  t8_hilbert_type_t type;
} t8_hilbert_t;

#endif /* T8_DQUAD_H */
