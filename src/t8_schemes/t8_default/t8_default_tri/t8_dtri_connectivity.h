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

/** \file t8_dtri_connectivity.h
 */

#ifndef T8_DTRI_CONNECTIVITY_H
#define T8_DTRI_CONNECTIVITY_H

#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_eclass.h>

T8_EXTERN_C_BEGIN ();

/** The spatial dimension */
#define T8_DTRI_DIM (2)

extern const int t8_tri_lut_type_vertex_to_cubevertex[2][4];
extern const int t8_tri_lut_cubevertex_to_num_adj[4];
extern const int t8_tri_lut_cubevertex_adj_to_type[4][2];
extern const int t8_tri_lut_cubevertex_adj_to_elementvertex[4][2];

/** Store the type of parent for each (cube-id,type) combination. */
extern const int t8_dtri_cid_type_to_parenttype[4][2];

/** Store the type of child for each (type,child number) combination,
  * where child number is the number in Bey order. */
extern const int t8_dtri_type_of_child[2][4];

/** Store the type of child for each (type,child number) combination,
  * where child number is the number in Morton order. */
extern const int t8_dtri_type_of_child_morton[2][4];

/** Store the Bey child number for each (Parent type,Morton child number) combination. */
extern const int t8_dtri_index_to_bey_number[2][4];

/** The anchor node of a child of a triangle T is the convex combination
 *  of T's anchor node x_0 and another node x_i of T.
 *  This array gives the index i in dependence of the Bey
 *  child id. */
extern const int t8_dtri_beyid_to_vertex[4];

/** Store the Bey child number for each (type,cube-id) combination. */
extern const int t8_dtri_type_cid_to_beyid[2][4];

/** Store the local index for each (parenttype,Bey child number) combination. */
extern const int t8_dtri_parenttype_beyid_to_Iloc[2][4];

/** Store the local index for each (type,cube-id) combination.*/
extern const int t8_dtri_type_cid_to_Iloc[2][4];

/** Store the type for each (parenttype,local Index) combination. */
extern const int t8_dtri_parenttype_Iloc_to_type[2][4];

/** Store the cube-id for each (parenttype,local Index) combination. */
extern const int t8_dtri_parenttype_Iloc_to_cid[2][4];

/** Store the indices of the corner of each face of a triangle. */
#define t8_dtri_face_corner t8_face_vertex_to_tree_vertex[T8_ECLASS_TRIANGLE]

/** Store the indices of the faces of each corner of a triangle. */
extern const int t8_dtri_corner_face[3][2];

T8_EXTERN_C_END ();

#endif /* T8_DTRI_CONNECTIVITY_H */
