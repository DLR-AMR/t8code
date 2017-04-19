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

/** \file t8_dtet_connectivity.h
 * TODO: come back later to see if it's worth having this separate file.
 */

#ifndef T8_DTET_CONNECTIVITY_H
#define T8_DTET_CONNECTIVITY_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** The spatial dimension */
#define T8_DTET_DIM (3)

/** Store the type of parent for each (cube-id,type) combination. */
extern const int    t8_dtet_cid_type_to_parenttype[8][6];

/** Store the type of child for each (type,child number) combination,
  * where child number is the number in Bey order. */
extern const int    t8_dtet_type_of_child[6][8];

/** Store the type of child for each (type,child number) combination,
  * where child number is the number in Morton order. */
extern const int    t8_dtet_type_of_child_morton[6][8];

/** Store the Bey child number for each (Parent type,Morton child number) combination. */
extern const int    t8_dtet_index_to_bey_number[6][8];

/** The anchor node of a child of a tetrahedron T is the convex combination
 *  of T's anchor node x_0 and another node x_i of T.
 *  This array gives the index i in dependence of the Bey
 *  child id. */
extern const int    t8_dtet_beyid_to_vertex[8];

/** Store the Bey child number for each (type,cube-id) combination. */
extern const int    t8_dtet_type_cid_to_beyid[6][8];

/** Store the local index for each (parenttype,Bey child number) combination. */
extern const int    t8_dtet_parenttype_beyid_to_Iloc[6][8];

/** Store the local index for each (type,cube-id) combination. */
extern const int    t8_dtet_type_cid_to_Iloc[6][8];

/** Store the type for each (parenttype,local Index) combination. */
extern const int    t8_dtet_parenttype_Iloc_to_type[6][8];

/** Store the cube-id for each (parenttype,local Index) combination. */
extern const int    t8_dtet_parenttype_Iloc_to_cid[6][8];

/** Store for each (type, face_index) the combination (category, type)
 *  of the respective boundary triangle.
 * I.e. {2, 1} means the boundary triangle is of category 2 and type 1.
 * The category determines how the coordinates of the triangle are computed
 * from the parent. \see t8_default_scheme_tet_c::t8_element_boundary
 */
extern const int    t8_dtet_type_face_to_boundary[6][4][2];

/** Store for each (type, face_index) the child_ids of the children of a tet of
 * the given type that share the given face.
 * I.e. {1,3} lists the child_ids of the children of a type 1 tetrahedron that have
 * a subface of face 3 of this tetrahedron.
 */
extern const int    t8_dtet_face_child_id_by_type[6][4][4];

T8_EXTERN_C_END ();

#endif /* T8_DTET_CONNECTIVITY_H */
