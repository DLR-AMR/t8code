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

/** \file t8_dquad_bits.h
 * Definitions of quad-specific functions.
 */

#ifndef T8_HILBERT_CONNECTIVITY_H
#define T8_HILBERT_CONNECTIVITY_H
#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Store the type of parent for each (cube-id,type) combination. */
extern const int8_t t8_hilbert_cubeid_type_to_parenttype[4][4];
extern const int8_t t8_hilbert_type_Iloc_to_childcubeid[4][4];
extern const int8_t t8_hilbert_type_Iloc_to_childtype[4][4];
extern const int8_t t8_hilbert_type_cubeid_to_Iloc[4][4];

T8_EXTERN_C_END ();

#endif /* T8_HILBERT_CONNECTIVITY_H */
