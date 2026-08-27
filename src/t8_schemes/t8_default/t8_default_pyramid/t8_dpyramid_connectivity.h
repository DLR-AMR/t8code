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

/** \file t8_dpyramid_connectivity.h
 * Definitions regarding the connectivity of Pyramids.
 */

#ifndef T8_DPYRAMID_CONNECTIVITY_H
#define T8_DPYRAMID_CONNECTIVITY_H

#include <t8.h>

/**The type of a pyramid depending on the parent pyramid and its local index
 *type = A(parent_type, local_index)
 */
extern const int t8_dpyramid_parenttype_Iloc_to_type[8][10];

/** The cube Id of a pyramid depending on its parenttype and local index
 * cube_id = A(parent_type, local_index)
 */
extern const int t8_dpyramid_parenttype_Iloc_to_cid[8][10];

/** The local ID of an element in a pyramid.
 * WARNING: The local ID is different, if the element is in a tet
 * Iloc = A(type, cube_id)*/
extern const int t8_dpyramid_type_cid_to_Iloc[8][8];

/** The type of the parent from of the cube id and its own type
 * parent_type = A(cube_id, type);
 */
extern const int t8_dpyramid_cid_type_to_parenttype[8][8];

/** The parenttype of a pyramid, computed by its own type and local ID
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * parenttype = A(type, local_index)
*/
extern const int t8_dpyramid_type_Iloc_to_parenttype[2][10];

/** The type of the parent of a pyramid, computed by its own type and cube-id
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * type = A(cube_id, parent_type)
*/
extern const int t8_dpyramid_type_cid_to_parenttype[2][8];

/** The number of local pyramid-siblings with lower id. This is computed with
 * help of the type of the parent and its own local id. A tetrahedron has
 * no pyramid-children, therefore this makes sense only for pyramidparents
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * num_siblings_with_lower_id = (parent_type, local_id)*/
extern const int t8_dpyramid_parenttype_iloc_pyra_w_lower_id[2][10];

/** Returns the face number of the neighbour touching the current pyramid.
 * The facenumber depends on the type of the pyramid
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * neigh_face_number = A(type, face)*/
extern const int t8_dpyramid_type_face_to_nface[2][5];

/** TODO: Documentation*/
extern const int t8_dpyramid_face_childid_to_is_inside[4][8];

/** The child ids of children touching a given face of a pyramid
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * child_id = A(type, face_num, child_id_at_face)
*/
extern const int t8_dpyramid_type_face_to_children_at_face[2][5][4];

/** Return the face-number of a children at a face of a pyramid
 * WARNING: The pyramid types 6 and 7 are encoded as 0 and 1! Can not be used for tets
 * face_number_of_child = A(type, face, child_id_at_face)
*/
extern const int t8_dpyramid_type_face_to_child_face[2][5][4];

/** Return the type of a boundary element which has a pyramid parent, depending on
 * the type of the boundary triangle and the face number of root
 * WARNING: The input type has to be the type of triangle
 * boundary_type = A(tri_type, face_number)*/
extern const int t8_dpyramid_tritype_rootface_to_pyratype[2][4];

/** Return the type of a boundary element which has a tet-parent, depending on
 * the type of the boundary triangle and the face number of root
 * WARNING: The input type has to be the type of triangle
 * tettype = A(tri_type, face_number_of_root)*/
extern const int t8_dpyramid_tritype_rootface_to_tettype[2][4];

/** Return the face number of a boundary element which has a tet-parent, depending
 * on the type of the boundary triangle and the the facenumber of root
 * WARNING: The input type has to be the type of triangle
 * face_number = A(tri_type, face_number_of_root)*/
extern const int t8_dpyramid_tritype_rootface_to_face[2][4];

/** Return the parent-type of a tetrahedron that has a pyramid as a parent
 * using the type of the tetrahedron and the cube_id of the tetrahedron
 * WARNING: The input type has to be the type of a tetrahedron
 * parent_type = A(tet_type, cube_id)
*/
extern const int t8_dtet_type_cid_to_pyramid_parenttype[6][8];

/** The corner numbers of the face of a pyramid
 * corner_number = A(pyramid_type, face, corner_number_at_face)
*/
extern const int t8_dpyramid_face_corner[2][5][4];

#endif  // T8_DPYRAMID_CONNECTIVITY_H
