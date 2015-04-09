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

/** \file t8_dtet_connectivity.h
 */

#ifndef T8_DTET_CONNECTIVITY_H
#define T8_DTET_CONNECTIVITY_H

/** The spatial dimension */
#define T8_DTET_DIM 3
/** The number of faces of a tetrahedron */
#define T8_DTET_FACES (T8_DTET_DIM +1)
/** The number of children of a tetrahedron
 *
 * also the nmber of corners */
#define T8_DTET_CHILDREN 8

/** Store the type of parent for each (cube-id,type) combination. */
extern const int                 t8_dtet_cid_type_to_parenttype[8][6];

/** Store the type of child for each (type,child number) combination. */
extern const int                 t8_dtet_type_of_child[6][8];

#endif /* T8_DTET_CONNECTIVITY_H */
