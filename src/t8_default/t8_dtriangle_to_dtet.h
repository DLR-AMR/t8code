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

/** \file t8_dtet_to_dtriangle.h
 */

#ifndef T8_DTRIANGLE_TO_DTET_H
#define T8_DTRIANGLE_TO_DTET_H

#define T8_DTR_TO_DTET

/* redefine macros */
#define T8_DTRIANGLE_MAXLEVEL T8_DTET_MAXLEVEL
#define T8_DTRIANGLE_ROOT_LEN T8_DTET_ROOT_LEN
#define T8_DTRIANGLE_DIM T8_DTET_DIM
#define T8_DTRIANGLE_CHILDREN T8_DTET_CHILDREN

/* redefine types */
#define t8_dtriangle_coord_t t8_dtet_coord_t
#define t8_dtriangle_type_t t8_dtet_type_t
#define t8_dtriangle_t       t8_dtet_t
#define t8_dtriangle_cube_id_t t8_dtet_cube_id_t

/* external variables */
#define t8_dtriangle_cid_type_to_parenttype t8_dtet_cid_type_to_parenttype
#define t8_dtriangle_type_of_child t8_dtet_type_of_child

/* functions in d8_dtriangle_bits.h */
#define t8_dtriangle_is_equal t8_dtet_is_equal
#define t8_dtriangle_parent t8_dtet_parent
#define t8_dtriangle_compute_coords t8_dtet_compute_coords
#define t8_dtriangle_child t8_dtet_child
#define t8_dtriangle_sibling t8_dtet_sibling
#define t8_dtriangle_face_neighbour t8_dtet_face_neighbour
#define t8_dtriangle_is_sibling t8_dtet_is_sibling
#define t8_dtriangle_is_parent t8_dtet_is_parent
#define t8_dtriangle_is_ancestor t8_dtet_is_ancestor

#endif /* T8_DTET_TO_DTRIANGLE_H */
