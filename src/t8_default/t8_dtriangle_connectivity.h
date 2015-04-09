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

/** \file t8_dtriangle_connectivity.h
 */

#ifndef T8_DTRIANGLE_CONNECTIVITY_H
#define T8_DTRIANGLE_CONNECTIVITY_H

#include "t8_dtriangle.h"

T8_EXTERN_C_BEGIN ();

/** The spatial dimension */
#define T8_DTRIANGLE_DIM 2

/** Store the type of parent for each (cube-id,type) combination. */
extern const int    t8_dtriangle_cid_type_to_parenttype[4][2];

/** Store the type of child for each (type,child number) combination. */
extern const int    t8_dtriangle_type_of_child[2][4];

T8_EXTERN_C_END ();

#endif /* T8_DTRIANGLE_CONNECTIVITY_H */
