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

/** \file t8_cmesh_trees.h
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_PART_TREE_H
#define T8_CMESH_PART_TREE_H

#include <t8.h>
#include <t8_cmesh.h>
#include "t8_cmesh_types.h"

T8_EXTERN_C_BEGIN ();

/* Interface for the data layout of the corse trees.
 *
 * The layout is the same for replicated and partitioned meshes.
 * Each process stores a meta array of data arrays. In the replicated case this meta
 * array has only one entry wheras in the partitioned case there is one data array for
 * each processor from which local trees were received in the last partition step
 * (and only one meta array if the cmesh arised from a partitioned commit).
 *
 * Each dara arrays stores the local trees, the ghosts, face neighbor information
 * of the ghosts, face neihbor information of the trees and the attributes of the trees.
 * Furthermore we store for each tree and for each ghost to which data array they belong to.
 * So the data looks like:
 *
 * M_0:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes |
 * M_1:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes |
 *  .         .        .          .            .               .
 * M_n:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes |
 *
 * tree_to_proc:  | 0 | 0 | 1 | ... | n |  these are just random examples here
 * ghost_to_proc: | 0 | 1 | 2 | ... | n |
 *
 *
 * Each tree T stores an offset to its Tree faces, such that (char*)&T + offset is
 * a pointer to the faces array.
 * The same holds for the ghost.
 * Also each tree stores the number of attributes and an offset relative to itself
 * to the first attribute entry of that tree.
 *
 * Tree faces:
 *
 * The data of Tree faces looks for each tree:
 *
 * | Treeid1 Treeid2  ... | ttf1 ttf2 ... | padding |
 *
 * Where padding is a number of unused bytes that makes the whole block a multiple
 * of 4 Bytes.
 * Treeid is a t8_locidx_t storing the local tree id for local tree neighbors and
 * the local ghost id + num_local_trees for ghost neighbors.
 * For the encoding of ttf (tree to face) see \ref t8_ctree_struct_t, ttf entries are int8_t
 * and the offset of ttf1 can be calculated from the Tree faces offset and the
 * class of the tree.
 *
 *  TODO: Should we really keep the padding?
 *
 * Ghost faces:
 *
 * | Treeid1 Treeid2 ... | ttf1 ttf2 ... | padding |
 *
 * Where padding is a number of unused bytes that makes the whole block a multiple
 * of 4 Bytes.
 * Treeidj is a t8_gloidx_t storing the global tree id for all neighbors.
 * For the encoding of ttf (tree to face) see \ref t8_ctree_struct_t, ttf entries are int8_t
 * and the offset of ttf1 can be calculated from the Tree faces offset and the
 * class of the tree.
 * Tree attributes:
 *
 * The data of Tree attributes looks like:
 *
 * TODO: This description seems incomplete. Why is it not Attij_data?
 *
 * | Att00_descr | Att01_descr | ... | Att10_descr | ... | Attrend_descr | Att1_data | Att2_data | ... |
 *                TODO: maybe insert padding here ||
 * Where Attij_descr is a descriptor of the j-th attribute data of tree i storing
 * - an offset to Atti_data starting from Atti0_descr
 * - package id of the attribute (int)
 * - key of the attribute (int)
 * The data type is t8_attribute_info_struct_t
 *
 * Attrend_descr only stores the offset of the end of this attributes block
 * (like an imaginary very last attribute);
 * using this info the size of each attribute can be computed as the difference
 * of the sizes of two consecutive attributes.
 *
 * padding is a number of nonused bytes to make the size of the descr block
 * a multiple of four.
 *
 *  TODO: maybe padding after the last Att_data is useful too
 *
 */

/* allocate a t8_cmesh_tree struct and allocate memory for its entries.
 * No memory for ctrees or ghosts is allocated here */
/* TODO: document */

/* Given a tree return the beginning of its attributes block */
#define T8_TREE_FIRST_ATT(t) ((char *)(t) + (t)->att_offset)
/* Given a tree and an index i return the i-th attribute index of that tree */
#define T8_TREE_ATTR_INFO(t,i) ((t8_attribute_info_struct_t *) \
  ((char*)(t) + (t)->att_offset + \
  (i) * sizeof (t8_attribute_info_struct_t)))
/* Given a tree and an attribute info return the attribute */
#define T8_TREE_ATTR(t,ai) (T8_TREE_FIRST_ATT(t) + (ai)->attribute_offset)
/* Given a tree return its face_neighbor array */
#define T8_TREE_FACE(t) ((char *) (t) + (t)->neigh_offset)
/* Given a tree return irs tree_to_face array */
#define T8_TREE_TTF(t) (T8_TREE_FACE(t) + \
  t8_eclass_num_faces[(t)->eclass] * sizeof(t8_locidx_t))
/* Given a ghost return its face_neighbor array */
#define T8_GHOST_FACE(g) T8_TREE_FACE(g)
/* Given a ghost return its tree_to_face array */
#define T8_GHOST_TTF(g) (int8_t *) (T8_GHOST_FACE(g) + \
  t8_eclass_num_faces[(g)->eclass] * sizeof(t8_gloidx_t))

void                t8_cmesh_trees_init (t8_cmesh_trees_t * ptrees,
                                         int num_procs, t8_locidx_t num_trees,
                                         t8_locidx_t num_ghosts);

#if 0
void                t8_cmesh_trees_init_part (t8_cmesh_trees_t trees,
                                              int proc,
                                              t8_locidx_t first_tree,
                                              t8_locidx_t last_tree,
                                              t8_locidx_t num_ghosts);
#endif

/* TODO: document */
t8_part_tree_t      t8_cmesh_trees_get_part (t8_cmesh_trees_t trees,
                                             int proc);

/* allocate the first_tree array of a given tree_part in a tree struct
 * with a given number of bytes */
/* !!! This does only allocate memory for the trees and ghosts
 *     not yet for the face data and the attributes. See below !!!
 */
void                t8_cmesh_trees_start_part (t8_cmesh_trees_t trees,
                                               int proc,
                                               t8_locidx_t first_tree,
                                               t8_locidx_t num_trees,
                                               t8_locidx_t first_ghost,
                                               t8_locidx_t num_ghosts);

/* TODO: document */
/* After all classes of trees and ghosts have been set and after the
 * number of tree attributes  was set and their total size (per tree)
 * stored temporarily in the att_offset variable
 * we grow the part array by the needed amount of memory and set the
 * offsets appropiately */
/* The workflow can be: call start_part, set tree and ghost classes maually, call
 * init_attributes, call finish_part, successively call add_attributes
 * and also set all face neighbors (TODO: write function)*/
void                t8_cmesh_trees_finish_part (t8_cmesh_trees_t trees,
                                                int proc);

/** Add a tree to a trees structure.
 * \param [in,out]  trees The trees structure to be updated.
 * \param [in]      tree_id The local id of the tree to be inserted.
 * \param [in]      proc  The mpirank of the process from which the tree was
 *                        received.
 * \param [in]      eclass The tree's element class.
 */
void                t8_cmesh_trees_add_tree (t8_cmesh_trees_t trees,
                                             t8_locidx_t ltree_id, int proc,
                                             t8_eclass_t eclass);

/** Add a ghost to a trees structure.
 * \param [in,out]  trees The trees structure to be updated.
 * \param [in]      ghost_index The local id of the ghost to be inserted.
 * \param [in]      tree_id The global index of the ghost.
 * \param [in]      proc  The mpirank of the process from which the ghost was
 *                        received.
 * \param [in]      eclass The ghost's element class.
 */
void                t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t ghost_index,
                                              t8_gloidx_t gtree_id, int proc,
                                              t8_eclass_t eclass);

/** Set all neighbor fields of all local trees and ghosts to boundary.
 * \param [in,out]  cmesh, The associated cmesh.
 * \param [in,out]  trees, The trees structure.
 */
void                t8_cmesh_trees_set_all_boundary (t8_cmesh_t cmesh,
                                                     t8_cmesh_trees_t trees);

/* TODO: This function return NULL if the tree is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific tree in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      tree  The local id of the tree.
 * \return                A pointer to the tree with local id \a tree.
 */
t8_ctree_t          t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees,
                                             t8_locidx_t tree);

/** Return a pointer to a specific tree in a trees struct plus pointers to
 * its face_neighbor and tree_to_face arrays.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      tree_id  The local id of the tree.
 * \param [out]     face_neigh If not NULL a pointer to the trees face_neighbor
 *                             array is stored here on return.
 * \param [out]     ttf        If not NULL a pointer to the trees tree_to_face
 *                             array is stored here on return.
 * \return                   A pointer to the tree with local id \a tree.
 */
t8_ctree_t          t8_cmesh_trees_get_tree_ext (t8_cmesh_trees_t trees,
                                                 t8_locidx_t tree_id,
                                                 t8_locidx_t ** face_neigh,
                                                 int8_t ** ttf);

/* TODO: This function return NULL if the ghost is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific ghost in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      ghost The local id of the ghost.
 * \return                A pointer to the ghost with local id \a ghost.
 */
t8_cghost_t         t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t ghost);

/** Return a pointer to a specific ghost in a trees struct plus pointers to
 * its face_neighbor and tree_to_face arrays.
 * \param [in]      trees The tress structure where the ghost is to be looked up.
 * \param [in]      ghost_id  The local id of the ghost.
 * \param [out]     face_neigh If not NULL a pointer to the ghosts face_neighbor
 *                             array is stored here on return.
 * \param [out]     ttf        If not NULL a pointer to the ghosts tree_to_face
 *                             array is stored here on return.
 * \return                   A pointer to the tree with local id \a tree.
 */
t8_cghost_t         t8_cmesh_trees_get_ghost_ext (t8_cmesh_trees_t trees,
                                                  t8_locidx_t ghost_id,
                                                  t8_gloidx_t ** face_neigh,
                                                  int8_t ** ttf);

/* TODO: document */
/* attr_bytes is the total size of all attributes of that tree */
void                t8_cmesh_trees_init_attributes (t8_cmesh_trees_t trees,
                                                    t8_locidx_t tree_id,
                                                    size_t num_attributes,
                                                    size_t attr_bytes);

/* TODO: document
 * sorts for each tree its attribute info objects, such that looking up
 * attributes is in O(log(A)) with A the number of attributes of that tree.
 * However, with this method we do not know the size of an attribute any longer,
 * this is something the user has to take care of */
void                t8_cmesh_trees_attribute_info_sort (t8_cmesh_trees_t
                                                        trees);

/* TODO: These need to be rewritten with package_id and key */
void               *t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees,
                                                  t8_topidx_t tree_id,
                                                  int package_id, int key);

/* TODO: document. total size of attributes of given tree */
size_t              t8_cmesh_trees_attribute_size (t8_ctree_t tree);

/* TODO: These need to be rewritten with package_id and key */
/* TODO: this uses char * and cmesh_set_attribute uses void *. Unify! */
/* attr_tree_index is index of attr in tree's attribute array.
 * We assume that the attributes are already sorted! */
void                t8_cmesh_trees_add_attribute (t8_cmesh_trees_t trees,
                                                  int proc,
                                                  t8_stash_attribute_struct_t
                                                  * attr, t8_locidx_t tree_id,
                                                  size_t index);

/* Print the trees,ghosts and their neighbors in ASCII.
 * This function is used for debugging purposes */
void                t8_cmesh_trees_print (t8_cmesh_t cmesh,
                                          t8_cmesh_trees_t trees);

int                 t8_cmesh_trees_is_equal (t8_cmesh_t cmesh,
                                             t8_cmesh_trees_t trees_a,
                                             t8_cmesh_trees_t trees_b);

void                t8_cmesh_trees_destroy (t8_cmesh_trees_t * trees);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PART_TREE_H */
