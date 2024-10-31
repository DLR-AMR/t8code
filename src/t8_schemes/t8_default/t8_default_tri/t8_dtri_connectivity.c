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

#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

/* clang-format off */
const int t8_dtri_cid_type_to_parenttype[4][2] = { 
  { 0, 1 }, 
  { 0, 0 }, 
  { 1, 1 }, 
  { 0, 1 } };

/* In dependence of a type x give the type of
 * the child with Bey number y */
const int t8_dtri_type_of_child[2][4] = { 
  { 0, 0, 0, 1 }, 
  { 1, 1, 1, 0 } };

/* in dependence of a type x give the type of
 * the child with Morton number y */
const int t8_dtri_type_of_child_morton[2][4] = { 
  { 0, 0, 1, 0 }, 
  { 1, 0, 1, 1 } };

/* Line b, row I gives the Bey child-id of
 * a Tet with Parent type b and local morton index I */
const int t8_dtri_index_to_bey_number[2][4] = { 
  { 0, 1, 3, 2 }, 
  { 0, 3, 1, 2 } };

const int t8_dtri_beyid_to_vertex[4] = { 0, 1, 2, 1 };

/* TODO: We us the next two tables after each other.
 *       We should replace this operation by a new table
 *          type + cid -> Iloc
 *       thus removing the Beyid from the operation.
 */

/* Line b, row c gives the Bey child-id of
 * a Tet with type b and cubeid c */
const int t8_dtri_type_cid_to_beyid[2][4] = { 
  { 0, 1, 3, 2 }, 
  { 0, 3, 1, 2 } };

/* Line b, row id gives the local index of
 * a Tet with type b and Bey child number id */
const int t8_dtri_parenttype_beyid_to_Iloc[2][4] = { 
  { 0, 1, 3, 2 }, 
  { 0, 2, 3, 1 } };

const int t8_dtri_type_cid_to_Iloc[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };

const int t8_dtri_parenttype_Iloc_to_type[2][4] = { 
  { 0, 0, 1, 0 }, 
  { 1, 0, 1, 1 } };

const int t8_dtri_parenttype_Iloc_to_cid[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };

const int t8_dtri_corner_face[3][2] = { 
  { 1, 2 }, 
  { 0, 2 }, 
  { 0, 1 } };

/* clang-format on*/
