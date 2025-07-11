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

/** \file t8_dvertex.hxx
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

/** The maximum refinement level allowed for a vertex. */
#define T8_DVERTEX_MAXLEVEL 255

template <>
struct t8_default_element <T8_ECLASS_VERTEX>
{
  uint8_t level;
};

typedef t8_default_element<T8_ECLASS_VERTEX> t8_dvertex_t;

#endif /* T8_DVERTEX_H */
