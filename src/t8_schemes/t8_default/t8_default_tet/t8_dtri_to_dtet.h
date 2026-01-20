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

/** \file t8_dtri_to_dtet.h 
 * In this file, we redefine macros, types,... of the code for triangles
 * to use the  default implementation also for tetrahedra in 3D to 
 * avoid code duplication.
 */

#ifndef T8_DTRI_TO_DTET_H
#define T8_DTRI_TO_DTET_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Define macro used in the default triangle implementation to check 
 * if we are in the case where the tri implementation is used for tets. 
 */
#define T8_DTRI_TO_DTET

/** Redefine macros. */
#define T8_DTRI_MAXLEVEL T8_DTET_MAXLEVEL           /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_ROOT_LEN T8_DTET_ROOT_LEN           /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_LEN T8_DTET_LEN                     /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_FACES T8_DTET_FACES                 /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_DIM T8_DTET_DIM                     /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_CHILDREN T8_DTET_CHILDREN           /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_FACES T8_DTET_FACES                 /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_FACE_CHILDREN T8_DTET_FACE_CHILDREN /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_CORNERS T8_DTET_CORNERS             /**< Wrapper of tri macro to tet.*/
#define T8_DTRI_NUM_TYPES T8_DTET_NUM_TYPES         /**< Wrapper of tri macro to tet.*/

/** Redefine types. */
#define t8_dtri_coord_t t8_dtet_coord_t     /**< Wrapper of tri type to tet.*/
#define t8_dtri_type_t t8_dtet_type_t       /**< Wrapper of tri type to tet.*/
#define t8_dtri_t t8_dtet_t                 /**< Wrapper of tri type to tet.*/
#define t8_dtri_cube_id_t t8_dtet_cube_id_t /**< Wrapper of tri type to tet.*/

/** External variables. */
#define t8_dtri_cid_type_to_parenttype t8_dtet_cid_type_to_parenttype   /**< Wrapper of tri variable to tet.*/
#define t8_dtri_type_of_child t8_dtet_type_of_child                     /**< Wrapper of tri variable to tet.*/
#define t8_dtri_type_of_child_morton t8_dtet_type_of_child_morton       /**< Wrapper of tri variable to tet.*/
#define t8_dtri_index_to_bey_number t8_dtet_index_to_bey_number         /**< Wrapper of tri variable to tet.*/
#define t8_dtri_beyid_to_vertex t8_dtet_beyid_to_vertex                 /**< Wrapper of tri variable to tet.*/
#define t8_dtri_type_cid_to_beyid t8_dtet_type_cid_to_beyid             /**< Wrapper of tri variable to tet.*/
#define t8_dtri_type_beyid_to_Iloc t8_dtet_type_beyid_to_Iloc           /**< Wrapper of tri variable to tet.*/
#define t8_dtri_parenttype_cid_to_Iloc t8_dtet_parenttype_cid_to_Iloc   /**< Wrapper of tri variable to tet.*/
#define t8_dtri_parenttype_Iloc_to_type t8_dtet_parenttype_Iloc_to_type /**< Wrapper of tri variable to tet.*/
#define t8_dtri_parenttype_Iloc_to_cid t8_dtet_parenttype_Iloc_to_cid   /**< Wrapper of tri variable to tet.*/
#define t8_dtri_type_cid_to_Iloc t8_dtet_type_cid_to_Iloc               /**< Wrapper of tri variable to tet.*/
#define t8_dtri_face_corner t8_dtet_face_corner                         /**< Wrapper of tri variable to tet.*/

/** Functions in t8_dtri_bits.h. */
#define t8_dtri_is_equal t8_dtet_is_equal                                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_copy t8_dtet_copy                                           /**< Wrapper of tri function to tet.*/
#define t8_dtri_compare t8_dtet_compare                                     /**< Wrapper of tri function to tet.*/
#define t8_dtri_equal t8_dtet_equal                                         /**< Wrapper of tri function to tet.*/
#define t8_dtri_parent t8_dtet_parent                                       /**< Wrapper of tri function to tet.*/
#define t8_dtri_ancestor t8_dtet_ancestor                                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_compute_all_coords t8_dtet_compute_all_coords               /**< Wrapper of tri function to tet.*/
#define t8_dtri_compute_integer_coords t8_dtet_compute_integer_coords       /**< Wrapper of tri function to tet.*/
#define t8_dtri_compute_vertex_ref_coords t8_dtet_compute_vertex_ref_coords /**< Wrapper of tri function to tet.*/
#define t8_dtri_compute_reference_coords t8_dtet_compute_reference_coords   /**< Wrapper of tri function to tet.*/
#define t8_dtri_child t8_dtet_child                                         /**< Wrapper of tri function to tet.*/
#define t8_dtri_childrenpv t8_dtet_childrenpv                               /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_familypv t8_dtet_is_familypv                             /**< Wrapper of tri function to tet.*/
#define t8_dtri_sibling t8_dtet_sibling                                     /**< Wrapper of tri function to tet.*/
#define t8_dtri_face_neighbour t8_dtet_face_neighbour                       /**< Wrapper of tri function to tet.*/
#define t8_dtri_nearest_common_ancestor t8_dtet_nearest_common_ancestor     /**< Wrapper of tri function to tet.*/
#define t8_dtri_children_at_face t8_dtet_children_at_face                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_face_child_face t8_dtet_face_child_face                     /**< Wrapper of tri function to tet.*/
#define t8_dtri_face_parent_face t8_dtet_face_parent_face                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_tree_face t8_dtet_tree_face                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_root_face_to_face t8_dtet_root_face_to_face                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_inside_root t8_dtet_is_inside_root                       /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_root_boundary t8_dtet_is_root_boundary                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_sibling t8_dtet_is_sibling                               /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_parent t8_dtet_is_parent                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_ancestor t8_dtet_is_ancestor                             /**< Wrapper of tri function to tet.*/
#define t8_dtri_linear_id t8_dtet_linear_id                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_linear_id_corner_desc t8_dtet_linear_id_corner_desc         /**< Wrapper of tri function to tet.*/
#define t8_dtri_init_linear_id t8_dtet_init_linear_id                       /**< Wrapper of tri function to tet.*/
#define t8_dtri_init_root t8_dtet_init_root                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_successor t8_dtet_successor                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_first_descendant t8_dtet_first_descendant                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_last_descendant t8_dtet_last_descendant                     /**< Wrapper of tri function to tet.*/
#define t8_dtri_corner_descendant t8_dtet_corner_descendant                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_predecessor t8_dtet_predecessor                             /**< Wrapper of tri function to tet.*/
#define t8_dtri_ancestor_id t8_dtet_ancestor_id                             /**< Wrapper of tri function to tet.*/
#define t8_dtri_child_id t8_dtet_child_id                                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_get_level t8_dtet_get_level                                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_is_valid t8_dtet_is_valid                                   /**< Wrapper of tri function to tet.*/
#define t8_dtri_init t8_dtet_init                                           /**< Wrapper of tri function to tet.*/
#define t8_dtri_init_linear_id_with_level t8_dtet_init_linear_id_with_level /**< Wrapper of tri function to tet.*/
#define t8_dtri_linear_id_with_level t8_dtet_linear_id_with_level           /**< Wrapper of tri function to tet.*/
#define t8_dtri_debug_print t8_dtet_debug_print                             /**< Wrapper of tri function to tet.*/
#define t8_dtri_element_pack t8_dtet_element_pack                           /**< Wrapper of tri function to tet.*/
#define t8_dtri_element_pack_size t8_dtet_element_pack_size                 /**< Wrapper of tri function to tet.*/
#define t8_dtri_element_unpack t8_dtet_element_unpack                       /**< Wrapper of tri function to tet.*/

T8_EXTERN_C_END ();

#endif /* T8_DTET_TO_DTRI_H */
