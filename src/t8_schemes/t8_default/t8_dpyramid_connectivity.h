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

/*The type of a pyramid depending on the parent pyramid and its local index
 *type = (parent_type, local_index)
 */
extern const int t8_dpyramid_parenttype_Iloc_to_type[8][10];

/*The cube Id of a pyramid depending on its parenttype and local index*/
extern const int t8_dpyramid_parenttype_Iloc_to_cid[8][10];

/* The local ID of an element in a pyramid. This is important!
 * The local ID is different, if the element is in a tet*/
extern const int           t8_dpyramid_type_cid_to_Iloc[8][8];

/* The type of the parent, dependant of the cube id and its own type
 */
extern const int t8_dpyramid_cid_type_to_parenttype[8][8];

/* The parenttype of a pyramid, computed by its won type and local ID*/
extern const int t8_dpyramid_type_Iloc_to_parenttype[2][10];

/*The type of the parent of a pyramid, computed by its own type and cube-id*/
extern const int t8_dpyramid_type_cid_to_parenttype[2][8];

/*The number of local pyramid-siblings with lower id. This is computed with
 * help of the type of the parent and its own local id. A tetrahedron has
 * no pyramid-children, therefore this makes sense only for pyramidparents*/
extern const int           t8_dpyramid_parenttype_iloc_pyra_w_lower_id[2][10];

/* Returns the facenumber of the neighbour touching the current pyramid.
 * The facenumber depends on the type of the pyramid*/
extern const int t8_dpyramid_type_face_to_nface[2][5];

#endif // T8_DPYRAMID_CONNECTIVITY_H

