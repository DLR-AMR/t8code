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

#ifndef T8_SIERPINSKI_CONNECTIVITY_H
#define T8_SIERPINSKI_CONNECTIVITY_H
#include <t8.h>

T8_EXTERN_C_BEGIN ();

extern const int8_t t8_sierpinski_parentimplicittype_cubeid_type_to_parenttype[2][4][2];
extern const int8_t t8_sierpinski_parentimplicittype_cubeid_type_to_Iloc[2][4][2];
extern const int8_t t8_sierpinski_parentimplicittype_parenttype_Iloc_to_childcubeid[2][2][4];
extern const int8_t t8_sierpinski_parentimplicittype_parenttype_Iloc_to_childtype[2][2][4];
extern const int8_t t8_sierpinski_ownimplicittype_type_vertex_to_cubevertex[2][2][3];

T8_EXTERN_C_END ();

#endif /* T8_SIERPINSKI_CONNECTIVITY_H */
