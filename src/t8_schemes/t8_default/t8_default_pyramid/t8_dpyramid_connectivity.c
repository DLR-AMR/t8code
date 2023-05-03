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

#include "t8_dpyramid_connectivity.h"
#include "t8_dpyramid.h"

/*The type of a child pyramid depending on the parent pyramid and its local index
 *type = (parent_type, local_index)
 */
const int           t8_dpyramid_type_Iloc_to_childtype[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CHILDREN] = {
  {0, 0, 2, 0, 1, 0, 1, 2, 3, 0},
  {1, 1, 3, 0, 1, 1, -1, -1, -1, -1},
  {2, 2, 3, 0, 2, 2, -1, -1, -1, -1},
  {3, 0, 1, 2, 3, 2, 3, 1, 3, 3}
};

/*The cube Id of a child pyramid depending on its parenttype and local index*/
const int
            t8_dpyramid_type_Iloc_to_childcubeid[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CHILDREN] = {
  {0, 1, 1, 2, 2, 3, 3, 3, 3, 7},
  {0, 1, 1, 5, 5, 7, -1, -1, -1, -1},
  {0, 2, 2, 6, 6, 7, -1, -1, -1, -1},
  {0, 4, 4, 4, 4, 5, 5, 6, 6, 7}
};

/* The local index of the element, dependant of the cube id and its own type
 */
const int           t8_dpyramid_type_cubeid_to_Iloc[T8_DPYRAMID_NUM_TYPES][1
                                                                           <<
                                                                           T8_DPYRAMID_DIM]
  = {
  {0, 1, 3, 5, 1, 3, 3, 9},
  {0, 1, 4, 6, 2, 4, 7, 5},
  {0, 2, 1, 7, 3, 5, 4, 5},
  {0, 2, 2, 8, 4, 6, 8, 9}
};

/* The type of the parent, dependant of the cube id and its own type
 */
const int
   
   
   
   
   
  t8_dpyramid_type_cubeid_to_parenttype[T8_DPYRAMID_NUM_TYPES][1 <<
                                                               T8_DPYRAMID_DIM]
  = {
  {0, 0, 0, 0, 3, 1, 2, 0},
  {1, 1, 0, 0, 3, 1, 3, 1},
  {2, 0, 2, 0, 3, 3, 2, 2},
  {3, 1, 2, 0, 3, 3, 3, 3}
};

const int
            t8_dpyramid_type_edge_equations[T8_DPYRAMID_NUM_EQUATIONS][2] = {
  {1, 2},
  {0, 2}
};

const int
            t8_dpyramid_type_vertex_dim_to_binary[T8_DPYRAMID_NUM_TYPES]
  [T8_DPYRAMID_MAX_CORNERS][T8_DPYRAMID_DIM] = {
  {
   {0, 0, 0},
   {1, 0, 0},
   {0, 1, 0},
   {1, 1, 0},
   {1, 1, 1}
   },
  {
   {0, 0, 0},
   {1, 0, 0},
   {1, 0, 1},
   {1, 1, 1},
   {-1, -1, -1}
   },
  {
   {0, 0, 0},
   {0, 1, 0},
   {0, 1, 1},
   {1, 1, 1},
   {-1, -1, -1}
   },
  {
   {0, 0, 0},
   {0, 0, 1},
   {1, 0, 1},
   {0, 1, 1},
   {1, 1, 1}
   }
};
