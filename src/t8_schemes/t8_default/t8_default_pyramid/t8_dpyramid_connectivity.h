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

#ifndef T8_DPYRAMID_CONNECTIVITY_H
#define T8_DPYRAMID_CONNECTIVITY_H

#include <t8.h>
#include "t8_dpyramid.h"

/* these two together enable computing the child*/

/**The type of a pyramid depending on the parent pyramid and its local index
 * child_type = A(type, local_index)
 */
extern const int   
  t8_dpyramid_type_Iloc_to_childtype[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CHILDREN];

/** The cube Id of a pyramid depending on its parenttype and local index
 * cube_id = A(type, local_index)
 */
extern const int   
  t8_dpyramid_type_Iloc_to_childcubeid[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CHILDREN];

/* shortcut to get own iloc*/

/** The local ID of an element in a pyramid.
 * Iloc = A(type, cube_id)*/
extern const int    t8_dpyramid_type_cubeid_to_Iloc[T8_DPYRAMID_NUM_TYPES][1
                                                                           <<
                                                                           T8_DPYRAMID_DIM];

/** The type of the parent of a pyramid, computed by its own type and cube-id
 * parent_type = A(type, cube_id)
*/
extern const int   
  t8_dpyramid_type_cubeid_to_parenttype[T8_DPYRAMID_NUM_TYPES][1 <<
                                                               T8_DPYRAMID_DIM];

extern const int   
  t8_dpyramid_type_edge_equations[T8_DPYRAMID_NUM_EQUATIONS][2];

extern const int   
  t8_dpyramid_type_vertex_dim_to_binary[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CORNERS][T8_DPYRAMID_DIM];

/**
extern const int    t8_dpyramid_type_cubeid_siblingid_to_siblingtype[T8_DPYRAMID_NUM_TYPES][1<<T8_DPYRAMID_DIM][T8_DPYRAMID_MAX_CHILDREN];
extern const int    t8_dpyramid_type_cubeid_siblingid_to_siblingcubeid[T8_DPYRAMID_NUM_TYPES][1<<T8_DPYRAMID_DIM][T8_DPYRAMID_MAX_CHILDREN];
*/

#endif // T8_DPYRAMID_CONNECTIVITY_H
