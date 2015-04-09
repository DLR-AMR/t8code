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

#ifndef T8_DTET_H
#define T8_DTET_H

#include <t8.h>

/** \file t8_dtet.h
 * TODO: document this.
 */

T8_EXTERN_C_BEGIN ();

/** The number of children that a tetrahedron is refined into. */
#define T8_DTET_CHILDREN 8

/** The maximum refinement level allowed for a tetrahedron. */
#define T8_DTET_MAXLEVEL 30

/** The length of the root tetrahedron in integer coordinates. */
#define T8_DTET_ROOT_LEN (1 << (T8_DTET_MAXLEVEL))

/** The length of a tetrahedron at a given level in integer coordinates. */
#define T8_DTET_LEN(l) (1 << (T8_DTET_MAXLEVEL - (l)))

/** The type of a tetrahedron designates its position relative to the surrounding cube. */
typedef int8_t      t8_dtet_type_t;

/** The coordinetes of a tetrahedron are integers relative to the maximum refinement. */
typedef int32_t     t8_dtet_coord_t;

/** This data type stores a tetrahedron. */
typedef struct t8_dtet
{
  /** We store the element class for compatibility with the pyramid. */
  int8_t              eclass;

  /** The refinement level of the tetrahedron relative to the root at level 0. */
  int8_t              level;

  /** Type of the tetrahedron in 0, ..., 5. */
  t8_dtet_type_t      type;

  t8_dtet_coord_t     x;        /**< The x integer coordinate of the anchor node. */
  t8_dtet_coord_t     y;        /**< The y integer coordinate of the anchor node. */
  t8_dtet_coord_t     z;        /**< The z integer coordinate of the anchor node. */
}
t8_dtet_t;

T8_EXTERN_C_END ();

#endif /* T8_DTET_H */
