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

/** \file t8_default_tet_bits.h
 */

#ifndef T8_DEFAULT_TET_BITS_H
#define T8_DEFAULT_TET_BITS_H

#include <t8_element.h>
#include "t8_default_tet.h"

void                t8_default_tet_parent (const t8_element_t * elem,
                                           t8_element_t * parent);

void                t8_default_tet_compute_coords (const t8_tet_t * t,
                                                   t8_tcoord_t
                                                   coordinates[4][3]);

void                t8_default_tet_child (const t8_element_t * elem,
                                          int childid, t8_element_t * child);

void                t8_default_tet_sibling (const t8_element_t * elem,
                                            int sibid,
                                            t8_element_t * sibling);

int                 t8_default_tet_face_neighbour (const t8_tet_t * t,
                                                   t8_tet_t * n, int face);

#endif /* T8_DEFAULT_TET_BITS_H */
