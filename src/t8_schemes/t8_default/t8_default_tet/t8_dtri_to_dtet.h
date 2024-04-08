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

/** \file t8_dtet_to_dtri.h
 */

#ifndef T8_DTRI_TO_DTET_H
#define T8_DTRI_TO_DTET_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

#define T8_DTRI_TO_DTET

/* redefine macros */
#define T8_DTRI_MAXLEVEL T8_DTET_MAXLEVEL
#define T8_DTRI_ROOT_LEN T8_DTET_ROOT_LEN
#define T8_DTRI_LEN T8_DTET_LEN
#define T8_DTRI_FACES T8_DTET_FACES
#define T8_DTRI_DIM T8_DTET_DIM
#define T8_DTRI_CHILDREN T8_DTET_CHILDREN
#define T8_DTRI_FACES T8_DTET_FACES
#define T8_DTRI_FACE_CHILDREN T8_DTET_FACE_CHILDREN
#define T8_DTRI_CORNERS T8_DTET_CORNERS
#define T8_DTRI_NUM_TYPES T8_DTET_NUM_TYPES

/* redefine types */
#define t8_dtri_coord_t t8_dtet_coord_t
#define t8_dtri_type_t t8_dtet_type_t
#define t8_dtri_t t8_dtet_t
#define t8_dtri_cube_id_t t8_dtet_cube_id_t

/* external variables */
#define t8_dtri_cid_type_to_parenttype t8_dtet_cid_type_to_parenttype
#define t8_dtri_type_of_child t8_dtet_type_of_child
#define t8_dtri_type_of_child_morton t8_dtet_type_of_child_morton
#define t8_dtri_index_to_bey_number t8_dtet_index_to_bey_number
#define t8_dtri_beyid_to_vertex t8_dtet_beyid_to_vertex
#define t8_dtri_type_cid_to_beyid t8_dtet_type_cid_to_beyid
#define t8_dtri_type_beyid_to_Iloc t8_dtet_type_beyid_to_Iloc
#define t8_dtri_parenttype_cid_to_Iloc t8_dtet_parenttype_cid_to_Iloc
#define t8_dtri_parenttype_Iloc_to_type t8_dtet_parenttype_Iloc_to_type
#define t8_dtri_parenttype_Iloc_to_cid t8_dtet_parenttype_Iloc_to_cid
#define t8_dtri_type_cid_to_Iloc t8_dtet_type_cid_to_Iloc
#define t8_dtri_face_corner t8_dtet_face_corner

/* functions in d8_dtri_bits.h */
#define t8_dtri_is_equal t8_dtet_is_equal
#define t8_dtri_copy t8_dtet_copy
#define t8_dtri_compare t8_dtet_compare
#define t8_dtri_equal t8_dtet_equal
#define t8_dtri_parent t8_dtet_parent
#define t8_dtri_ancestor t8_dtet_ancestor
#define t8_dtri_compute_all_coords t8_dtet_compute_all_coords
#define t8_dtri_compute_integer_coords t8_dtet_compute_integer_coords
#define t8_dtri_compute_vertex_ref_coords t8_dtet_compute_vertex_ref_coords
#define t8_dtri_compute_reference_coords t8_dtet_compute_reference_coords
#define t8_dtri_child t8_dtet_child
#define t8_dtri_childrenpv t8_dtet_childrenpv
#define t8_dtri_is_familypv t8_dtet_is_familypv
#define t8_dtri_sibling t8_dtet_sibling
#define t8_dtri_face_neighbour t8_dtet_face_neighbour
#define t8_dtri_nearest_common_ancestor t8_dtet_nearest_common_ancestor
#define t8_dtri_children_at_face t8_dtet_children_at_face
#define t8_dtri_face_child_face t8_dtet_face_child_face
#define t8_dtri_face_parent_face t8_dtet_face_parent_face
#define t8_dtri_tree_face t8_dtet_tree_face
#define t8_dtri_root_face_to_face t8_dtet_root_face_to_face
#define t8_dtri_is_inside_root t8_dtet_is_inside_root
#define t8_dtri_is_root_boundary t8_dtet_is_root_boundary
#define t8_dtri_is_sibling t8_dtet_is_sibling
#define t8_dtri_is_parent t8_dtet_is_parent
#define t8_dtri_is_ancestor t8_dtet_is_ancestor
#define t8_dtri_linear_id t8_dtet_linear_id
#define t8_dtri_linear_id_corner_desc t8_dtet_linear_id_corner_desc
#define t8_dtri_init_linear_id t8_dtet_init_linear_id
#define t8_dtri_init_root t8_dtet_init_root
#define t8_dtri_successor t8_dtet_successor
#define t8_dtri_first_descendant t8_dtet_first_descendant
#define t8_dtri_last_descendant t8_dtet_last_descendant
#define t8_dtri_corner_descendant t8_dtet_corner_descendant
#define t8_dtri_predecessor t8_dtet_predecessor
#define t8_dtri_ancestor_id t8_dtet_ancestor_id
#define t8_dtri_child_id t8_dtet_child_id
#define t8_dtri_get_level t8_dtet_get_level
#define t8_dtri_is_valid t8_dtet_is_valid
#define t8_dtri_init t8_dtet_init
#define t8_dtri_init_linear_id_with_level t8_dtet_init_linear_id_with_level
#define t8_dtri_linear_id_with_level t8_dtet_linear_id_with_level
#define t8_dtri_debug_print t8_dtet_debug_print
#define t8_dtri_element_pack t8_dtet_element_pack
#define t8_dtri_element_pack_size t8_dtet_element_pack_size
#define t8_dtri_element_unpack t8_dtet_element_unpack

T8_EXTERN_C_END ();

#endif /* T8_DTET_TO_DTRI_H */
