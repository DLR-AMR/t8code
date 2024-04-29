/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** We do not process this file in Doxygen since it overwrites types.
 *
 * Transform 2D \ref p4est routines into 3D \ref p8est routines.  This file can
 * be included from a .c file that has been written for 2D to turn it into a 3D
 * code with minor modifications.  We #define P4_TO_P8, which allows to compile
 * a source file twice, once for 2D without and once for 3D with including this
 * header, adding extra code for 3D if necessary.
 */

#ifndef P4EST_TO_P8EST_H
#define P4EST_TO_P8EST_H

#ifdef P4EST_H
#error "The include files p4est.h and p4est_to_p8est.h cannot be combined"
#endif
#define P4_TO_P8

#include <p4est_base.h>

/* redefine macros */
#define P4EST_ONDISK_FORMAT             P8EST_ONDISK_FORMAT
#define P4EST_DIM                       P8EST_DIM
#define P4EST_DIM_POW                   P8EST_DIM_POW
#define P4EST_FACES                     P8EST_FACES
#define P4EST_CHILDREN                  P8EST_CHILDREN
#define P4EST_HALF                      P8EST_HALF
#define P4EST_FTRANSFORM                P8EST_FTRANSFORM
#define P4EST_INSUL                     P8EST_INSUL
#define P4EST_ONLY_P8_LAND              P8EST_ONLY_P8_LAND
#define P4EST_ONLY_P8_COMMA             P8EST_ONLY_P8_COMMA
#define P4EST_STRING                    P8EST_STRING
#define P4EST_MAXLEVEL                  P8EST_MAXLEVEL
#define P4EST_QMAXLEVEL                 P8EST_QMAXLEVEL
#define P4EST_OLD_MAXLEVEL              P8EST_OLD_MAXLEVEL
#define P4EST_OLD_QMAXLEVEL             P8EST_OLD_QMAXLEVEL
#define P4EST_ROOT_LEN                  P8EST_ROOT_LEN
#define P4EST_QUADRANT_LEN              P8EST_QUADRANT_LEN
#define P4EST_QUADRANT_MASK             P8EST_QUADRANT_MASK
#define P4EST_LAST_OFFSET               P8EST_LAST_OFFSET
#define P4EST_QUADRANT_INIT             P8EST_QUADRANT_INIT
#define P4EST_LEAF_IS_FIRST_IN_TREE     P8EST_LEAF_IS_FIRST_IN_TREE

#ifdef P4EST_ENABLE_FILE_DEPRECATED

#define P4EST_FILE_MAGIC_NUMBER         P8EST_FILE_MAGIC_NUMBER
#define P4EST_FILE_METADATA_BYTES       P8EST_FILE_METADATA_BYTES
#define P4EST_FILE_MAGIC_BYTES          P8EST_FILE_MAGIC_BYTES
#define P4EST_FILE_VERSION_STR_BYTES    P8EST_FILE_VERSION_STR_BYTES
#define P4EST_FILE_ARRAY_METADATA_BYTES P8EST_FILE_ARRAY_METADATA_BYTES
#define P4EST_FILE_ARRAY_METADATA_CHARS P8EST_FILE_ARRAY_METADATA_CHARS
#define P4EST_FILE_BYTE_DIV             P8EST_FILE_BYTE_DIV
#define P4EST_FILE_MAX_NUM_PAD_BYTES    P8EST_FILE_MAX_NUM_PAD_BYTES
#define P4EST_FILE_USER_STRING_BYTES    P8EST_FILE_USER_STRING_BYTES
#define P4EST_FILE_FIELD_HEADER_BYTES   P8EST_FILE_FIELD_HEADER_BYTES
#define P4EST_FILE_MAX_GLOBAL_QUAD      P8EST_FILE_MAX_GLOBAL_QUAD
#define P4EST_FILE_MAX_BLOCK_SIZE       P8EST_FILE_MAX_BLOCK_SIZE
#define P4EST_FILE_MAX_FIELD_ENTRY_SIZE P8EST_FILE_MAX_FIELD_ENTRY_SIZE

#endif /* P4EST_ENABLE_FILE_DEPRECATED */

/* redefine enums */
#define P4EST_CONNECT_SELF              P8EST_CONNECT_SELF
#define P4EST_CONNECT_FACE              P8EST_CONNECT_FACE
#define P4EST_CONNECT_ALMOST            P8EST_CONNECT_ALMOST
#define P4EST_CONNECT_CORNER            P8EST_CONNECT_CORNER
#define P4EST_CONNECT_FULL              P8EST_CONNECT_FULL
#define P4EST_CONN_ENCODE_NONE          P8EST_CONN_ENCODE_NONE
#define P4EST_TRANSFER_COMM_SRC         P8EST_TRANSFER_COMM_SRC
#define P4EST_TRANSFER_COMM_DEST        P8EST_TRANSFER_COMM_DEST
#define P4EST_TRANSFER_COMM_SRC_DUP     P8EST_TRANSFER_COMM_SRC_DUP
#define P4EST_TRANSFER_COMM_DEST_DUP    P8EST_TRANSFER_COMM_DEST_DUP
#define P4EST_TRANSFER_COMM_EXTERNAL    P8EST_TRANSFER_COMM_EXTERNAL
#define P4EST_WRAP_NONE                 P8EST_WRAP_NONE
#define P4EST_WRAP_REFINE               P8EST_WRAP_REFINE
#define P4EST_WRAP_COARSEN              P8EST_WRAP_COARSEN

#ifdef P4EST_ENABLE_FILE_DEPRECATED

#define P4EST_FILE_ERR_SUCCESS          P8EST_FILE_ERR_SUCCESS
#define P4EST_FILE_ERR_FILE             P8EST_FILE_ERR_FILE
#define P4EST_FILE_ERR_NOT_SAME         P8EST_FILE_ERR_NOT_SAME
#define P4EST_FILE_ERR_AMODE            P8EST_FILE_ERR_AMODE
#define P4EST_FILE_ERR_NO_SUCH_FILE     P8EST_FILE_ERR_NO_SUCH_FILE
#define P4EST_FILE_ERR_FILE_EXIST       P8EST_FILE_ERR_FILE_EXIST
#define P4EST_FILE_ERR_BAD_FILE         P8EST_FILE_ERR_BAD_FILE
#define P4EST_FILE_ERR_ACCESS           P8EST_FILE_ERR_ACCESS
#define P4EST_FILE_ERR_NO_SPACE         P8EST_FILE_ERR_NO_SPACE
#define P4EST_FILE_ERR_QUOTA            P8EST_FILE_ERR_QUOTA
#define P4EST_FILE_ERR_READ_ONLY        P8EST_FILE_ERR_READ_ONLY
#define P4EST_FILE_ERR_IN_USE           P8EST_FILE_ERR_IN_USE
#define P4EST_FILE_ERR_IO               P8EST_FILE_ERR_IO
#define P4EST_FILE_ERR_FORMAT           P8EST_FILE_ERR_FORMAT
#define P4EST_FILE_ERR_SECTION_TYPE     P8EST_FILE_ERR_SECTION_TYPE
#define P4EST_FILE_ERR_CONN             P8EST_FILE_ERR_CONN
#define P4EST_FILE_ERR_P4EST            P8EST_FILE_ERR_P8EST
#define P4EST_FILE_ERR_IN_DATA          P8EST_FILE_ERR_IN_DATA
#define P4EST_FILE_ERR_COUNT            P8EST_FILE_ERR_COUNT
#define P4EST_FILE_ERR_UNKNOWN          P8EST_FILE_ERR_UNKNOWN
#define P4EST_FILE_ERR_LASTCODE         P8EST_FILE_ERR_LASTCODE

#endif

/* redefine types */
#ifdef P4EST_BACKWARD_DEALII
#define p4est_balance_type_t            p8est_balance_type_t
#endif
#define p4est_connect_type_t            p8est_connect_type_t
#define p4est_connectivity_encode_t     p8est_connectivity_encode_t
#define p4est_connectivity_t            p8est_connectivity_t
#define p4est_corner_transform_t        p8est_corner_transform_t
#define p4est_corner_info_t             p8est_corner_info_t
#define p4est_neighbor_transform_t      p8est_neighbor_transform_t
#define p4est_geometry_t                p8est_geometry_t
#define p4est_geometry_destroy_t        p8est_geometry_destroy_t
#define p4est_t                         p8est_t
#define p4est_tree_t                    p8est_tree_t
#define p4est_quadrant_t                p8est_quadrant_t
#define p4est_inspect_t                 p8est_inspect_t
#define p4est_position_t                p8est_position_t
#define p4est_init_t                    p8est_init_t
#define p4est_refine_t                  p8est_refine_t
#define p4est_coarsen_t                 p8est_coarsen_t
#define p4est_weight_t                  p8est_weight_t
#define p4est_ghost_t                   p8est_ghost_t
#define p4est_ghost_exchange_t          p8est_ghost_exchange_t
#define p4est_indep_t                   p8est_indep_t
#define p4est_nodes_t                   p8est_nodes_t
#define p4est_lid_t                     p8est_lid_t
#define p4est_lnodes_t                  p8est_lnodes_t
#define p4est_lnodes_code_t             p8est_lnodes_code_t
#define p4est_lnodes_rank_t             p8est_lnodes_rank_t
#define p4est_lnodes_buffer_t           p8est_lnodes_buffer_t
#define p4est_iter_volume_t             p8est_iter_volume_t
#define p4est_iter_volume_info_t        p8est_iter_volume_info_t
#define p4est_iter_face_t               p8est_iter_face_t
#define p4est_iter_face_info_t          p8est_iter_face_info_t
#define p4est_iter_face_side_t          p8est_iter_face_side_t
#define p4est_iter_corner_t             p8est_iter_corner_t
#define p4est_iter_corner_side_t        p8est_iter_corner_side_t
#define p4est_iter_corner_info_t        p8est_iter_corner_info_t
#define p4est_mesh_params_t             p8est_mesh_params_t
#define p4est_search_query_t            p8est_search_query_t
#define p4est_search_local_t            p8est_search_local_t
#define p4est_search_reorder_t          p8est_search_reorder_t
#define p4est_search_partition_t        p8est_search_partition_t
#define p4est_search_all_t              p8est_search_all_t
#define p4est_build                     p8est_build
#define p4est_build_t                   p8est_build_t
#define p4est_transfer_comm_t           p8est_transfer_comm_t
#define p4est_transfer_context_t        p8est_transfer_context_t
#define p4est_mesh_t                    p8est_mesh_t
#define p4est_mesh_face_neighbor_t      p8est_mesh_face_neighbor_t
#define p4est_wrap_t                    p8est_wrap_t
#define p4est_wrap_leaf_t               p8est_wrap_leaf_t
#define p4est_wrap_flags_t              p8est_wrap_flags_t
#define p4est_wrap_params_t             p8est_wrap_params_t
#define p4est_vtk_context_t             p8est_vtk_context_t
#define p4est_file_context_t            p8est_file_context_t
#define p4est_file_section_metadata_t   p8est_file_section_metadata_t

/* redefine external variables */
#define p4est_face_corners              p8est_face_corners
#define p4est_face_dual                 p8est_face_dual
#define p4est_corner_faces              p8est_corner_faces
#define p4est_corner_face_corners       p8est_corner_face_corners
#define p4est_child_corner_faces        p8est_child_corner_faces
#define P4EST_DATA_UNINITIALIZED        P8EST_DATA_UNINITIALIZED

/* functions in p4est_connectivity */
#define p4est_connectivity_face_neighbor_face_corner    \
        p8est_connectivity_face_neighbor_face_corner
#define p4est_connectivity_face_neighbor_corner         \
        p8est_connectivity_face_neighbor_corner
#define p4est_connectivity_memory_used  p8est_connectivity_memory_used
#define p4est_connectivity_new          p8est_connectivity_new
#define p4est_connectivity_new_brick    p8est_connectivity_new_brick
#define p4est_connectivity_new_periodic p8est_connectivity_new_periodic
#define p4est_connectivity_new_twotrees p8est_connectivity_new_twotrees
#define p4est_connectivity_new_byname   p8est_connectivity_new_byname
#define p4est_connectivity_new_copy     p8est_connectivity_new_copy
#define p4est_connectivity_bcast        p8est_connectivity_bcast
#define p4est_connectivity_destroy      p8est_connectivity_destroy
#define p4est_connectivity_set_attr     p8est_connectivity_set_attr
#define p4est_connectivity_is_valid     p8est_connectivity_is_valid
#define p4est_connectivity_is_equal     p8est_connectivity_is_equal
#define p4est_connectivity_sink         p8est_connectivity_sink
#define p4est_connectivity_deflate      p8est_connectivity_deflate
#define p4est_connectivity_save         p8est_connectivity_save
#define p4est_connectivity_source       p8est_connectivity_source
#define p4est_connectivity_inflate      p8est_connectivity_inflate
#define p4est_connectivity_load         p8est_connectivity_load
#define p4est_connectivity_complete     p8est_connectivity_complete
#define p4est_connectivity_reduce       p8est_connectivity_reduce
#define p4est_expand_face_transform     p8est_expand_face_transform
#define p4est_find_face_transform       p8est_find_face_transform
#define p4est_find_corner_transform     p8est_find_corner_transform
#define p4est_neighbor_transform_coordinates \
        p8est_neighbor_transform_coordinates
#define p4est_neighbor_transform_coordinates_reverse \
        p8est_neighbor_transform_coordinates_reverse
#define p4est_connectivity_get_neighbor_transforms \
        p8est_connectivity_get_neighbor_transforms
#define p4est_corner_array_index        p8est_corner_array_index
#define p4est_connectivity_reorder      p8est_connectivity_reorder
#define p4est_connectivity_reorder_newid                \
        p8est_connectivity_reorder_newid
#define p4est_connectivity_permute      p8est_connectivity_permute
#define p4est_connectivity_join_faces   p8est_connectivity_join_faces
#define p4est_connectivity_is_equivalent p8est_connectivity_is_equivalent
#define p4est_connectivity_read_inp_stream p8est_connectivity_read_inp_stream
#define p4est_connectivity_read_inp     p8est_connectivity_read_inp

/* functions in p4est */
#define p4est_qcoord_to_vertex          p8est_qcoord_to_vertex
#define p4est_memory_used               p8est_memory_used
#define p4est_revision                  p8est_revision
#define p4est_new                       p8est_new
#define p4est_destroy                   p8est_destroy
#define p4est_copy                      p8est_copy
#define p4est_reset_data                p8est_reset_data
#define p4est_refine                    p8est_refine
#define p4est_coarsen                   p8est_coarsen
#define p4est_balance                   p8est_balance
#define p4est_partition                 p8est_partition
#define p4est_checksum                  p8est_checksum
#define p4est_checksum_partition        p8est_checksum_partition
#define p4est_save                      p8est_save
#define p4est_load                      p8est_load
#define p4est_connect_type_int          p8est_connect_type_int
#define p4est_connect_type_string       p8est_connect_type_string
#define p4est_tree_array_index          p8est_tree_array_index
#define p4est_quadrant_array_index      p8est_quadrant_array_index
#define p4est_quadrant_array_push_copy  p8est_quadrant_array_push_copy
#define p4est_quadrant_array_push       p8est_quadrant_array_push
#define p4est_quadrant_mempool_alloc    p8est_quadrant_mempool_alloc
#define p4est_quadrant_list_pop         p8est_quadrant_list_pop

/* functions in p4est_extended */
#define p4est_replace_t                 p8est_replace_t
#define p4est_lid_compare               p8est_lid_compare
#define p4est_lid_is_equal              p8est_lid_is_equal
#define p4est_lid_init                  p8est_lid_init
#define p4est_lid_set_zero              p8est_lid_set_zero
#define p4est_lid_set_one               p8est_lid_set_one
#define p4est_lid_set_uint64            p8est_lid_set_uint64
#define p4est_lid_chk_bit               p8est_lid_chk_bit
#define p4est_lid_set_bit               p8est_lid_set_bit
#define p4est_lid_copy                  p8est_lid_copy
#define p4est_lid_add                   p8est_lid_add
#define p4est_lid_sub                   p8est_lid_sub
#define p4est_lid_bitwise_neg           p8est_lid_bitwise_neg
#define p4est_lid_bitwise_or            p8est_lid_bitwise_or
#define p4est_lid_bitwise_and           p8est_lid_bitwise_and
#define p4est_lid_shift_right           p8est_lid_shift_right
#define p4est_lid_shift_left            p8est_lid_shift_left
#define p4est_lid_add_inplace           p8est_lid_add_inplace
#define p4est_lid_sub_inplace           p8est_lid_sub_inplace
#define p4est_lid_bitwise_or_inplace    p8est_lid_bitwise_or_inplace
#define p4est_lid_bitwise_and_inplace   p8est_lid_bitwise_and_inplace
#define p4est_quadrant_linear_id_ext128 p8est_quadrant_linear_id_ext128
#define p4est_quadrant_set_morton_ext128 p8est_quadrant_set_morton_ext128
#define p4est_new_ext                   p8est_new_ext
#define p4est_mesh_new_ext              p8est_mesh_new_ext
#define p4est_mesh_new_params           p8est_mesh_new_params
#define p4est_mesh_params_init          p8est_mesh_params_init
#define p4est_copy_ext                  p8est_copy_ext
#define p4est_refine_ext                p8est_refine_ext
#define p4est_coarsen_ext               p8est_coarsen_ext
#define p4est_balance_ext               p8est_balance_ext
#define p4est_balance_subtree_ext       p8est_balance_subtree_ext
#define p4est_partition_ext             p8est_partition_ext
#define p4est_partition_for_coarsening  p8est_partition_for_coarsening
#define p4est_save_ext                  p8est_save_ext
#define p4est_load_ext                  p8est_load_ext
#define p4est_source_ext                p8est_source_ext

#ifdef P4EST_ENABLE_FILE_DEPRECATED

#define p4est_file_open_read_ext        p8est_file_open_read_ext
#define p4est_file_read_field_ext       p8est_file_read_field_ext

#endif

/* functions in p4est_iterate */
#define p4est_iterate                   p8est_iterate
#define p4est_iterate_ext               p8est_iterate_ext
#define p4est_iter_fside_array_index    p8est_iter_fside_array_index
#define p4est_iter_fside_array_index_int p8est_iter_fside_array_index_int
#define p4est_iter_cside_array_index    p8est_iter_cside_array_index
#define p4est_iter_cside_array_index_int p8est_iter_cside_array_index_int

/* functions in p4est_points */
#define p4est_new_points                p8est_new_points

/* functions in p4est_bits */
#define p4est_quadrant_pad              p8est_quadrant_pad
#define p4est_quadrant_print            p8est_quadrant_print
#define p4est_quadrant_is_equal         p8est_quadrant_is_equal
#define p4est_quadrant_overlaps         p8est_quadrant_overlaps
#define p4est_quadrant_is_equal_piggy   p8est_quadrant_is_equal_piggy
#define p4est_quadrant_compare          p8est_quadrant_compare
#define p4est_coordinates_compare       p8est_coordinates_compare
#define p4est_quadrant_disjoint         p8est_quadrant_disjoint
#define p4est_quadrant_compare_piggy    p8est_quadrant_compare_piggy
#define p4est_quadrant_compare_local_num p8est_quadrant_compare_local_num
#define p4est_quadrant_equal_fn         p8est_quadrant_equal_fn
#define p4est_quadrant_hash_fn          p8est_quadrant_hash_fn
#define p4est_node_equal_piggy_fn       p8est_node_equal_piggy_fn
#define p4est_node_hash_piggy_fn        p8est_node_hash_piggy_fn
#define p4est_node_clamp_inside         p8est_node_clamp_inside
#define p4est_node_unclamp              p8est_node_unclamp
#define p4est_node_to_quadrant          p8est_node_to_quadrant
#define p4est_quadrant_contains_node    p8est_quadrant_contains_node
#define p4est_quadrant_ancestor_id      p8est_quadrant_ancestor_id
#define p4est_quadrant_child_id         p8est_quadrant_child_id
#define p4est_coordinates_is_inside_root p8est_coordinates_is_inside_root
#define p4est_quadrant_is_inside_root   p8est_quadrant_is_inside_root
#define p4est_quadrant_is_inside_3x3    p8est_quadrant_is_inside_3x3
#define p4est_quadrant_is_outside_face  p8est_quadrant_is_outside_face
#define p4est_quadrant_is_outside_corner p8est_quadrant_is_outside_corner
#define p4est_quadrant_is_node          p8est_quadrant_is_node
#define p4est_coordinates_is_valid      p8est_coordinates_is_valid
#define p4est_quadrant_is_valid         p8est_quadrant_is_valid
#define p4est_quadrant_is_extended      p8est_quadrant_is_extended
#define p4est_quadrant_is_sibling       p8est_quadrant_is_sibling
#define p4est_quadrant_is_sibling_D     p8est_quadrant_is_sibling_D
#define p4est_quadrant_is_family        p8est_quadrant_is_family
#define p4est_quadrant_is_familyv       p8est_quadrant_is_familyv
#define p4est_quadrant_is_familypv      p8est_quadrant_is_familypv
#define p4est_quadrant_is_parent        p8est_quadrant_is_parent
#define p4est_quadrant_is_parent_D      p8est_quadrant_is_parent_D
#define p4est_quadrant_is_ancestor      p8est_quadrant_is_ancestor
#define p4est_quadrant_is_ancestor_D    p8est_quadrant_is_ancestor_D
#define p4est_quadrant_is_next          p8est_quadrant_is_next
#define p4est_quadrant_is_next_D        p8est_quadrant_is_next_D
#define p4est_quadrant_overlaps_tree    p8est_quadrant_overlaps_tree
#define p4est_quadrant_is_inside_tree   p8est_quadrant_is_inside_tree
#define p4est_quadrant_is_first_last    p8est_quadrant_is_first_last
#define p4est_quadrant_enlarge_first    p8est_quadrant_enlarge_first
#define p4est_quadrant_enlarge_last     p8est_quadrant_enlarge_last
#define p4est_quadrant_ancestor         p8est_quadrant_ancestor
#define p4est_quadrant_parent           p8est_quadrant_parent
#define p4est_quadrant_sibling          p8est_quadrant_sibling
#define p4est_quadrant_child            p8est_quadrant_child
#define p4est_quadrant_face_neighbor    p8est_quadrant_face_neighbor
#define p4est_quadrant_face_neighbor_extra p8est_quadrant_face_neighbor_extra
#define p4est_quadrant_half_face_neighbors p8est_quadrant_half_face_neighbors
#define p4est_quadrant_all_face_neighbors p8est_quadrant_all_face_neighbors
#define p4est_quadrant_corner_neighbor  p8est_quadrant_corner_neighbor
#define p4est_quadrant_corner_neighbor_extra    \
        p8est_quadrant_corner_neighbor_extra
#define p4est_quadrant_half_corner_neighbor     \
        p8est_quadrant_half_corner_neighbor
#define p4est_quadrant_corner_node      p8est_quadrant_corner_node
#define p4est_quadrant_children         p8est_quadrant_children
#define p4est_quadrant_childrenv        p8est_quadrant_childrenv
#define p4est_quadrant_childrenpv       p8est_quadrant_childrenpv
#define p4est_quadrant_first_descendant p8est_quadrant_first_descendant
#define p4est_quadrant_last_descendant  p8est_quadrant_last_descendant
#define p4est_quadrant_corner_descendant p8est_quadrant_corner_descendant
#define p4est_nearest_common_ancestor   p8est_nearest_common_ancestor
#define p4est_nearest_common_ancestor_D p8est_nearest_common_ancestor_D
#define p4est_coordinates_transform_face        \
        p8est_coordinates_transform_face
#define p4est_quadrant_transform_face   p8est_quadrant_transform_face
#define p4est_quadrant_touches_corner   p8est_quadrant_touches_corner
#define p4est_quadrant_transform_corner p8est_quadrant_transform_corner
#define p4est_quadrant_shift_corner     p8est_quadrant_shift_corner
#define p4est_quadrant_linear_id        p8est_quadrant_linear_id
#define p4est_quadrant_set_morton       p8est_quadrant_set_morton
#define p4est_quadrant_successor        p8est_quadrant_successor
#define p4est_quadrant_predecessor      p8est_quadrant_predecessor
#define p4est_quadrant_srand            p8est_quadrant_srand
#define p4est_neighbor_transform_quadrant       \
        p8est_neighbor_transform_quadrant
#define p4est_neighbor_transform_quadrant_reverse    \
        p8est_neighbor_transform_quadrant_reverse
#define p4est_quadrant_is_ancestor_face          \
        p8est_quadrant_is_ancestor_face
#define p4est_quadrant_is_ancestor_corner        \
        p8est_quadrant_is_ancestor_corner

/* functions in p4est_search */
#define p4est_find_partition            p8est_find_partition
#define p4est_find_lower_bound          p8est_find_lower_bound
#define p4est_find_higher_bound         p8est_find_higher_bound
#define p4est_find_quadrant_cumulative  p8est_find_quadrant_cumulative
#define p4est_split_array               p8est_split_array
#define p4est_find_range_boundaries     p8est_find_range_boundaries
#define p4est_search                    p8est_search
#define p4est_search_local              p8est_search_local
#define p4est_search_reorder            p8est_search_reorder
#define p4est_search_partition          p8est_search_partition
#define p4est_search_partition_gfx      p8est_search_partition_gfx
#define p4est_search_partition_gfp      p8est_search_partition_gfp
#define p4est_search_all                p8est_search_all
#define p4est_build_new                 p8est_build_new
#define p4est_build_init_add            p8est_build_init_add
#define p4est_build_add                 p8est_build_add
#define p4est_build_complete            p8est_build_complete

/* functions in p4est_algorithms */
#define p4est_quadrant_mempool_new      p8est_quadrant_mempool_new
#define p4est_quadrant_init_data        p8est_quadrant_init_data
#define p4est_quadrant_free_data        p8est_quadrant_free_data
#define p4est_quadrant_checksum         p8est_quadrant_checksum
#define p4est_quadrant_in_range         p8est_quadrant_in_range
#define p4est_tree_is_sorted            p8est_tree_is_sorted
#define p4est_tree_is_linear            p8est_tree_is_linear
#define p4est_tree_is_almost_sorted     p8est_tree_is_almost_sorted
#define p4est_tree_is_complete          p8est_tree_is_complete
#define p4est_tree_print                p8est_tree_print
#define p4est_is_equal                  p8est_is_equal
#define p4est_quadrant_copy             p8est_quadrant_copy
#define p4est_is_valid                  p8est_is_valid
#define p4est_tree_compute_overlap      p8est_tree_compute_overlap
#define p4est_tree_uniqify_overlap      p8est_tree_uniqify_overlap
#define p4est_tree_remove_nonowned      p8est_tree_remove_nonowned
#define p4est_complete_region           p8est_complete_region
#define p4est_complete_subtree          p8est_complete_subtree
#define p4est_balance_subtree           p8est_balance_subtree
#define p4est_balance_border            p8est_balance_border
#define p4est_linearize_tree            p8est_linearize_tree
#define p4est_next_nonempty_process     p8est_next_nonempty_process
#define p4est_partition_correction      p8est_partition_correction
#define p4est_partition_for_coarsening  p8est_partition_for_coarsening
#define p4est_partition_given           p8est_partition_given
#define p4est_quadrant_on_face_boundary p8est_quadrant_on_face_boundary

/* functions in p4est_communication */
#define p4est_comm_parallel_env_assign  p8est_comm_parallel_env_assign
#define p4est_comm_parallel_env_duplicate p8est_comm_parallel_env_duplicate
#define p4est_comm_parallel_env_release p8est_comm_parallel_env_release
#define p4est_comm_parallel_env_replace p8est_comm_parallel_env_replace
#define p4est_comm_parallel_env_get_info p8est_comm_parallel_env_get_info
#define p4est_comm_parallel_env_is_null p8est_comm_parallel_env_is_null
#define p4est_comm_parallel_env_reduce  p8est_comm_parallel_env_reduce
#define p4est_comm_parallel_env_reduce_ext p8est_comm_parallel_env_reduce_ext
#define p4est_comm_count_quadrants      p8est_comm_count_quadrants
#define p4est_comm_global_partition     p8est_comm_global_partition
#define p4est_comm_global_first_quadrant p8est_comm_global_first_quadrant
#define p4est_comm_count_pertree        p8est_comm_count_pertree
#define p4est_comm_is_empty             p8est_comm_is_empty
#define p4est_comm_is_empty_gfq         p8est_comm_is_empty_gfq
#define p4est_comm_is_empty_gfp         p8est_comm_is_empty_gfp
#define p4est_comm_is_contained         p8est_comm_is_contained
#define p4est_comm_is_owner             p8est_comm_is_owner
#define p4est_comm_is_owner_gfp         p8est_comm_is_owner_gfp
#define p4est_comm_find_owner           p8est_comm_find_owner
#define p4est_comm_tree_info            p8est_comm_tree_info
#define p4est_comm_neighborhood_owned   p8est_comm_neighborhood_owned
#define p4est_comm_sync_flag            p8est_comm_sync_flag
#define p4est_comm_checksum             p8est_comm_checksum
#define p4est_transfer_fixed            p8est_transfer_fixed
#define p4est_bsearch_partition         p8est_bsearch_partition
#define p4est_transfer_fixed_begin      p8est_transfer_fixed_begin
#define p4est_transfer_fixed_end        p8est_transfer_fixed_end
#define p4est_transfer_custom           p8est_transfer_custom
#define p4est_transfer_custom_begin     p8est_transfer_custom_begin
#define p4est_transfer_custom_end       p8est_transfer_custom_end
#define p4est_transfer_items            p8est_transfer_items
#define p4est_transfer_items_begin      p8est_transfer_items_begin
#define p4est_transfer_items_end        p8est_transfer_items_end
#define p4est_transfer_end              p8est_transfer_end

/* functions in p4est_io */
#define p4est_deflate_quadrants         p8est_deflate_quadrants
#define p4est_inflate                   p8est_inflate
#define p4est_inflate_null              p8est_inflate_null

#ifdef P4EST_ENABLE_FILE_DEPRECATED

#define p4est_file_open_create          p8est_file_open_create
#define p4est_file_open_append          p8est_file_open_append
#define p4est_file_open_read            p8est_file_open_read
#define p4est_file_write_block          p8est_file_write_block
#define p4est_file_read_block           p8est_file_read_block
#define p4est_file_write_field          p8est_file_write_field
#define p4est_file_read_field           p8est_file_read_field
#define p4est_file_info                 p8est_file_info
#define p4est_file_error_string         p8est_file_error_string
#define p4est_file_write_p4est          p8est_file_write_p8est
#define p4est_file_read_p4est           p8est_file_read_p8est
#define p4est_file_write_connectivity   p8est_file_write_connectivity
#define p4est_file_read_connectivity    p8est_file_read_connectivity
#define p4est_file_close                p8est_file_close

#endif /* P4EST_ENABLE_FILE_DEPRECATED */

/* functions in p4est_geometry */
#define p4est_geometry_destroy          p8est_geometry_destroy
#define p4est_geometry_new_connectivity p8est_geometry_new_connectivity
#define p4est_geometry_connectivity_X   p8est_geometry_connectivity_X

/* functions in p4est_vtk */
#define p4est_vtk_context_new           p8est_vtk_context_new
#define p4est_vtk_context_destroy       p8est_vtk_context_destroy
#define p4est_vtk_context_set_geom      p8est_vtk_context_set_geom
#define p4est_vtk_context_set_scale     p8est_vtk_context_set_scale
#define p4est_vtk_context_set_continuous p8est_vtk_context_set_continuous
#define p4est_vtk_write_file            p8est_vtk_write_file
#define p4est_vtk_write_header          p8est_vtk_write_header
#define p4est_vtk_write_header_ho       p8est_vtk_write_header_ho
#define p4est_vtk_write_cell_dataf      p8est_vtk_write_cell_dataf
#define p4est_vtk_write_cell_datav      p8est_vtk_write_cell_datav
#define p4est_vtk_write_cell_data       p8est_vtk_write_cell_data
#define p4est_vtk_write_point_dataf     p8est_vtk_write_point_dataf
#define p4est_vtk_write_point_data      p8est_vtk_write_point_data
#define p4est_vtk_write_footer          p8est_vtk_write_footer

/* functions in p4est_ghost */
#define p4est_quadrant_find_owner       p8est_quadrant_find_owner
#define p4est_ghost_memory_used         p8est_ghost_memory_used
#define p4est_ghost_new                 p8est_ghost_new
#define p4est_ghost_new_local           p8est_ghost_new_local
#define p4est_ghost_destroy             p8est_ghost_destroy
#define p4est_ghost_exchange_data       p8est_ghost_exchange_data
#define p4est_ghost_exchange_data_begin p8est_ghost_exchange_data_begin
#define p4est_ghost_exchange_data_end   p8est_ghost_exchange_data_end
#define p4est_ghost_exchange_custom     p8est_ghost_exchange_custom
#define p4est_ghost_exchange_custom_begin p8est_ghost_exchange_custom_begin
#define p4est_ghost_exchange_custom_end p8est_ghost_exchange_custom_end
#define p4est_ghost_exchange_custom_levels p8est_ghost_exchange_custom_levels
#define p4est_ghost_exchange_custom_levels_begin        \
        p8est_ghost_exchange_custom_levels_begin
#define p4est_ghost_exchange_custom_levels_end  \
        p8est_ghost_exchange_custom_levels_end
#define p4est_ghost_bsearch             p8est_ghost_bsearch
#define p4est_ghost_contains            p8est_ghost_contains
#define p4est_ghost_is_valid            p8est_ghost_is_valid
#define p4est_face_quadrant_exists      p8est_face_quadrant_exists
#define p4est_quadrant_exists           p8est_quadrant_exists
#define p4est_is_balanced               p8est_is_balanced
#define p4est_ghost_checksum            p8est_ghost_checksum
#define p4est_ghost_expand              p8est_ghost_expand

/* functions in p4est_nodes */
#define p4est_nodes_new                 p8est_nodes_new
#define p4est_nodes_destroy             p8est_nodes_destroy
#define p4est_nodes_is_valid            p8est_nodes_is_valid

/* functions in p4est_lnodes */
#define p4est_lnodes_new                p8est_lnodes_new
#define p4est_lnodes_destroy            p8est_lnodes_destroy
#define p4est_ghost_support_lnodes      p8est_ghost_support_lnodes
#define p4est_ghost_expand_by_lnodes    p8est_ghost_expand_by_lnodes
#define p4est_partition_lnodes          p8est_partition_lnodes
#define p4est_partition_lnodes_detailed p8est_partition_lnodes_detailed
#define p4est_lnodes_decode             p8est_lnodes_decode
#define p4est_lnodes_share_owned_begin  p8est_lnodes_share_owned_begin
#define p4est_lnodes_share_owned_end    p8est_lnodes_share_owned_end
#define p4est_lnodes_share_owned        p8est_lnodes_share_owned
#define p4est_lnodes_share_all_begin    p8est_lnodes_share_all_begin
#define p4est_lnodes_share_all_end      p8est_lnodes_share_all_end
#define p4est_lnodes_share_all          p8est_lnodes_share_all
#define p4est_lnodes_buffer_destroy     p8est_lnodes_buffer_destroy
#define p4est_lnodes_rank_array_index   p8est_lnodes_rank_array_index
#define p4est_lnodes_rank_array_index_int       \
        p8est_lnodes_rank_array_index_int
#define p4est_lnodes_global_index       p8est_lnodes_global_index

/* functions in p4est_mesh */
#define p4est_mesh_memory_used          p8est_mesh_memory_used
#define p4est_mesh_new                  p8est_mesh_new
#define p4est_mesh_destroy              p8est_mesh_destroy
#define p4est_mesh_get_quadrant         p8est_mesh_get_quadrant
#define p4est_mesh_get_neighbors        p8est_mesh_get_neighbors
#define p4est_mesh_quadrant_cumulative  p8est_mesh_quadrant_cumulative
#define p4est_mesh_face_neighbor_init   p8est_mesh_face_neighbor_init
#define p4est_mesh_face_neighbor_init2  p8est_mesh_face_neighbor_init2
#define p4est_mesh_face_neighbor_next   p8est_mesh_face_neighbor_next
#define p4est_mesh_face_neighbor_data   p8est_mesh_face_neighbor_data

/* functions in p4est_balance */
#define p4est_balance_seeds_face        p8est_balance_seeds_face
#define p4est_balance_seeds_corner      p8est_balance_seeds_corner
#define p4est_balance_seeds             p8est_balance_seeds

/* functions in p4est_wrap */
#define p4est_wrap_params_init          p8est_wrap_params_init
#define p4est_wrap_new_conn             p8est_wrap_new_conn
#define p4est_wrap_new_p4est            p8est_wrap_new_p8est
#define p4est_wrap_new_p4est_params     p8est_wrap_new_p8est_params
#define p4est_wrap_new_brick            p8est_wrap_new_brick
#define p4est_wrap_new_world            p8est_wrap_new_world
#define p4est_wrap_new_ext              p8est_wrap_new_ext
#define p4est_wrap_new_params           p8est_wrap_new_params
#define p4est_wrap_new_copy             p8est_wrap_new_copy
#define p4est_wrap_destroy              p8est_wrap_destroy
#define p4est_wrap_set_hollow           p8est_wrap_set_hollow
#define p4est_wrap_set_coarsen_delay    p8est_wrap_set_coarsen_delay
#define p4est_wrap_set_partitioning     p8est_wrap_set_partitioning
#define p4est_wrap_get_ghost            p8est_wrap_get_ghost
#define p4est_wrap_get_mesh             p8est_wrap_get_mesh
#define p4est_wrap_mark_refine          p8est_wrap_mark_refine
#define p4est_wrap_mark_coarsen         p8est_wrap_mark_coarsen
#define p4est_wrap_adapt                p8est_wrap_adapt
#define p4est_wrap_partition            p8est_wrap_partition
#define p4est_wrap_complete             p8est_wrap_complete
#define p4est_wrap_leaf_next            p8est_wrap_leaf_next
#define p4est_wrap_leaf_first           p8est_wrap_leaf_first

/* functions in p4est_plex */
#define p4est_get_plex_data             p8est_get_plex_data
#define p4est_get_plex_data_ext         p8est_get_plex_data_ext

/* functions in p4est_connrefine */
#define p4est_connectivity_refine       p8est_connectivity_refine

#endif /* !P4EST_TO_P8EST_H */
