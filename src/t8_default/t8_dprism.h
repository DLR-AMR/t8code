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

#ifndef T8_PRISM_H
#define T8_PRISM_H

/** \file t8_dprism.h
 * TODO: document this.
 */

#include <t8.h>
#include <./t8_default/t8_dline.h>
#include <./t8_default/t8_dtri.h>

T8_EXTERN_C_BEGIN ();

/** The number of children that a line is refined into. */
#define T8_DPRISM_CHILDREN 8

/** The number of faces of a line. */
#define T8_DPRISM_FACES 5

/** The maximum refinement level allowed for a line. */
#define T8_DPRISM_MAXLEVEL 30

/** The length of the root line in integer coordinates. */
#define T8_DPRISM_ROOT_LEN (1 << (T8_DPRISM_MAXLEVEL))

/** The length of a line at a given level in integer coordinates. */
#define T8_DPRISM_LEN(l) (1 << (T8_DPRISM_MAXLEVEL - (l)))

typedef int32_t     t8_dprism_coord_t;

typedef struct t8_dprism
{
  t8_dline_t          line;     /*z coordinate + level */
  t8_dtri_t           tri;      /*x,y coordinate + level */
}
t8_dprism_t;

T8_EXTERN_C_END ();

#endif /* T8_DPRISM_H */
