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

/* Interface for the data layout of the coarse trees.
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
 * TODO: would be more logical to switch Ghost and Tree faces
 *
 * M_0:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes | Ghost attributes
 * M_1:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes | Ghost attributes
 *  .         .        .          .            .               .
 * M_n:   | Trees | Ghosts | Ghost faces | Tree faces | Tree attributes | Ghost attributes
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
 * Ghost faces:
 *
 * | Treeid1 Treeid2 ... | ttf1 ttf2 ... | padding |
 *
 * Where padding is a number of unused bytes that makes the whole block a multiple
 * of sizeof (void*) Bytes.
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
 *  TODO: We may also need padding between the attributes.
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

/* Given a ghost return the beginning of its attribute block */
#define T8_GHOST_FIRST_ATT(g) T8_TREE_FIRST_ATT (g)

/* Given a ghost and an index i return the i-th attribute index of that ghost */
#define T8_GHOST_ATTR_INFO(g,i) T8_TREE_ATTR_INFO (g, i)

/* Given a ghost and an attribute info return the attribute */
#define T8_GHOST_ATTR(g,ai) T8_TREE_ATTR(g,ai)

/* Given a ghost return its face_neighbor array */
#define T8_GHOST_FACE(g) T8_TREE_FACE(g)

/* Given a ghost return its tree_to_face array */
#define T8_GHOST_TTF(g) (int8_t *) (T8_GHOST_FACE(g) + \
  t8_eclass_num_faces[(g)->eclass] * sizeof(t8_gloidx_t))

/** This struct is an entry of the trees global_id to local_id
 * hash table for ghost trees. */
typedef struct
{
  t8_gloidx_t         global_id;/**< The global id */
  t8_locidx_t         local_id; /**< The local id */
} t8_trees_glo_lo_hash_t;

/** Initialize a trees structure and allocate its parts.
 * This function allocates the from_procs array without filling it, it
 * also allocates the tree_to_proc and ghost_to_proc arrays.
 * No memory for trees or ghosts is allocated.
 * \param [in,ou      ptrees   The trees structure to be initialized.
 * \param [in]        num_procs The number of entries of its from_proc array
 *                              (can be different for each process).
 * \param [in]        num_trees The number of trees that will be stored in this
 *                              structure.
 * \param [in]        num_ghosts The number of ghosts that will be stored in this
 *                              structure.
 */
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

/** Return one part of a specified tree array.
 * \param [in]        trees   The tree array to be queried
 * \param [in]        proc    An index specifying the part to be returned.
 * \return                    The part number \a proc of \a trees.
 */
t8_part_tree_t      t8_cmesh_trees_get_part (t8_cmesh_trees_t trees,
                                             int proc);

/* !!! This does only allocate memory for the trees and ghosts
 *     not yet for the face data and the attributes. See below !!!
 */
/** Allocate the first_tree array of a given tree_part in a tree struct
 *  with a given number of trees and ghosts.
 *  This function allocates the memory for the trees and the ghosts
 *  but not for their face neighbor entries or attributes. These must
 *  be allocated later when the eclasses of the trees and ghosts are known
 *  \ref t8_cmesh_trees_finish_part.
 *  \param [in,out]         trees   The trees structure to be updated.
 *  \param [in]             proc    The index of the part to be updated.
 *  \param [in]             lfirst_tree The local id of the first tree of that part.
 *  \param [in]             num_trees The number of trees of that part.
 *  \param [in]             lfirst_ghost The local id of the first ghost of that part.
 *  \param [in]             num_ghosts The number of ghosts of that part.
 *  \param [in]             alloc   If true then the first_tree array is allocated for
 *                          the number of trees and ghosts.
 *                          When a cmesh is copied we do not want this, so in we pass alloc = 0 then.
 */
void                t8_cmesh_trees_start_part (t8_cmesh_trees_t trees,
                                               int proc,
                                               t8_locidx_t lfirst_tree,
                                               t8_locidx_t num_trees,
                                               t8_locidx_t lfirst_ghost,
                                               t8_locidx_t num_ghosts,
                                               int alloc);

/** After all classes of trees and ghosts have been set and after the
 * number of tree attributes  was set and their total size (per tree)
 * stored temporarily in the att_offset variable
 * we grow the part array by the needed amount of memory and set the
 * offsets appropiately.
 * The workflow should be: call \ref t8_cmesh_trees_start_part,
 * set tree and ghost classes maually via \ref t8_cmesh_trees_add_tree
 * and \ref t8_cmesh_trees_add_ghost, call
 * \ref t8_cmesh_trees_init_attributes, then call this function.
 * Afterwards successively call \ref t8_cmesh_trees_add_attribute for
 * each attribute and
 * also set all face neighbors (TODO: write function).
 * \param [in,out]        trees The trees structure to be updated.
 * \param [in]            proc  The number of the part to be finished.
 */
void                t8_cmesh_trees_finish_part (t8_cmesh_trees_t trees,
                                                int proc);

/** Copy the tree_to_proc and ghost_to_proc arrays of one tree structure to
 * another one.
 * \param [in,out]      trees_dest    The destination trees structure.
 * \param [in]          trees_src     The source trees structure.
 * \param [in]          lnum_trees    The total number of trees stored in \a trees_src.
 * \param [in]          lnum_ghosts    The total number of ghosts stored in \a trees_src.
 */
void                t8_cmesh_trees_copy_toproc (t8_cmesh_trees_t trees_dest,
                                                t8_cmesh_trees_t trees_src,
                                                t8_locidx_t lnum_trees,
                                                t8_locidx_t lnum_ghosts);

/** Copy the trees array from one part to another.
 * \param [in,out]      trees_dest    The trees struct of the destination part.
 * \param [in]          part_dest     The index of the destination part. Must be initialized
 *                                    by \ref t8_cmesh_trees_start_part with alloc = 0.
 * \param [in]          trees_src     The trees struct of the source part.
 * \param [in]          part_src      The index of the destination part.
 *                                    Must be a valid part, thus \ref t8_cmesh_trees_finish_part
 *                                    must have been called.
 */
void                t8_cmesh_trees_copy_part (t8_cmesh_trees_t trees_dest,
                                              int part_dest,
                                              t8_cmesh_trees_t trees_src,
                                              int part_src);

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
 * \param [in]      ghost_index The index in the part array of the ghost to be inserted.
 * \param [in]      tree_id The global index of the ghost.
 * \param [in]      proc  The mpirank of the process from which the ghost was
 *                        received.
 * \param [in]      eclass The ghost's element class.
 * \param [in]      num_local_trees The number of local trees in the cmesh.
 */
void                t8_cmesh_trees_add_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t lghost_index,
                                              t8_gloidx_t gtree_id, int proc,
                                              t8_eclass_t eclass,
                                              t8_locidx_t num_local_trees);

/** Set all neighbor fields of all local trees and ghosts to boundary.
 * \param [in,out]  cmesh, The associated cmesh.
 * \param [in,out]  trees, The trees structure.
 * A face f of tree t counts as boundary if the face-neighbor is also t
 * at face f.
 */
void                t8_cmesh_trees_set_all_boundary (t8_cmesh_t cmesh,
                                                     t8_cmesh_trees_t trees);

void                t8_cmesh_trees_get_part_data (t8_cmesh_trees_t trees,
                                                  int proc,
                                                  t8_locidx_t * first_tree,
                                                  t8_locidx_t * num_trees,
                                                  t8_locidx_t * first_ghost,
                                                  t8_locidx_t * num_ghosts);

/* TODO: This function returns NULL if the tree is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific tree in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      ltree  The local id of the tree.
 * \return                A pointer to the tree with local id \a tree.
 */
t8_ctree_t          t8_cmesh_trees_get_tree (t8_cmesh_trees_t trees,
                                             t8_locidx_t ltree);

/** Return a pointer to a specific tree in a trees struct plus pointers to
 * its face_neighbor and tree_to_face arrays.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      ltree_id  The local id of the tree.
 * \param [out]     face_neigh If not NULL a pointer to the trees face_neighbor
 *                             array is stored here on return.
 * \param [out]     ttf        If not NULL a pointer to the trees tree_to_face
 *                             array is stored here on return.
 * \return                   A pointer to the tree with local id \a tree.
 */
t8_ctree_t          t8_cmesh_trees_get_tree_ext (t8_cmesh_trees_t trees,
                                                 t8_locidx_t ltree_id,
                                                 t8_locidx_t ** face_neigh,
                                                 int8_t ** ttf);

/** Given a coarse tree and a face number, return the local id of the neighbor tree.
 * \param [in]      tree.     The coarse tree.
 * \param [in]      face.     The face number.
 * \return                    The local id of the neighbor tree. */
t8_locidx_t         t8_cmesh_trees_get_face_neighbor (t8_ctree_t tree,
                                                      int face);

/* TODO: This function returns NULL if the ghost is not present.
 *       So far no error checking is done here. */
/** Return a pointer to a specific ghost in a trees struct.
 * \param [in]      trees The tress structure where the tree is to be looked up.
 * \param [in]      lghost The local id of the ghost.
 * \return                A pointer to the ghost with local id \a ghost.
 */
t8_cghost_t         t8_cmesh_trees_get_ghost (t8_cmesh_trees_t trees,
                                              t8_locidx_t lghost);

/** Return a pointer to a specific ghost in a trees struct plus pointers to
 * its face_neighbor and tree_to_face arrays.
 * \param [in]      trees The trees structure where the ghost is to be looked up.
 * \param [in]      lghost_id  The local id of the ghost.
 * \param [out]     face_neigh If not NULL a pointer to the ghosts face_neighbor
 *                             array is stored here on return.
 * \param [out]     ttf        If not NULL a pointer to the ghosts tree_to_face
 *                             array is stored here on return.
 * \return                   A pointer to the tree with local id \a tree.
 */
t8_cghost_t         t8_cmesh_trees_get_ghost_ext (t8_cmesh_trees_t trees,
                                                  t8_locidx_t lghost_id,
                                                  t8_gloidx_t ** face_neigh,
                                                  int8_t ** ttf);

/** Given the global tree id of a ghost tree in a trees structure,
 * return its local ghost id.
 * \param [in]      trees   The trees structure.
 * \param [in]      global_id A global tree id.
 * \return                  The local id of the tree \a global_id if it is a ghost
 *                          in \a trees. A negative number if it isn't.
 *                          The local id is a number l with
 *                          num_local_trees <= \a l < num_local_trees + num_ghosts
 */
t8_locidx_t         t8_cmesh_trees_get_ghost_local_id (t8_cmesh_trees_t trees,
                                                       t8_gloidx_t global_id);

/* TODO: document.
 * returns the complete size in bytes needed to store all information */
size_t              t8_cmesh_trees_size (t8_cmesh_trees_t trees);

/** For one tree in a trees structure set the number of attributes
 *  and temporarily store the total size of all of this tree's attributes.
 *  This temporary value is used in \ref t8_cmesh_trees_finish_part.
 * \param [in,out]        trees The trees structure to be updated.
 * \param [in]            ltree_id The local id of one tree in \a trees.
 * \param [in]            num_attributes The number of attributes of this tree.
 * \param [in]            attr_bytes The total number of bytes of all attributes
 *                                   of this tree.
 */
void                t8_cmesh_trees_init_attributes (t8_cmesh_trees_t trees,
                                                    t8_locidx_t ltree_id,
                                                    size_t num_attributes,
                                                    size_t attr_bytes);

#if 0
/* TODO: document
 * TODO: Is this function needed?
 * sorts for each tree its attribute info objects, such that looking up
 * attributes is in O(log(A)) with A the number of attributes of that tree.
 * However, with this method we do not know the size of an attribute any longer,
 * this is something the user has to take care of */
void                t8_cmesh_trees_attribute_info_sort (t8_cmesh_trees_t
                                                        trees);
#endif

/** Return an attribute that is stored at a tree.
 *  \param [in]       trees   The trees structure.
 *  \param [in]       ltree_id  The local id of the tree whose attribute is querid.
 *  \param [in]       package_id The package identifier of the attribute.
 *  \param [in]       key       The key of the attribute within all attributes of
 *                              the same package identifier.
 *  \param [out]      size      If not NULL, the size (in bytes) of the attribute
 *                              will be stored here.
 *  \param [in]       is_ghost  If true, then \a ltree_id is interpreted as the local_id
 *                              of a ghost.
 *  \return           A pointer to the queried attribute, NULL if the attribute
 *                    does not exist.
 */
void               *t8_cmesh_trees_get_attribute (t8_cmesh_trees_t trees,
                                                  t8_locidx_t ltree_id,
                                                  int package_id, int key,
                                                  size_t * size,
                                                  int is_ghost);

/** Return the total size of all attributes stored at a specified tree.
 * \param [in]        tree  A tree structure.
 * \return            The total size (in bytes) of the attributes of \a tree.
 */
size_t              t8_cmesh_trees_attribute_size (t8_ctree_t tree);

/** Return the total size of all attributes stored at a specified ghost.
 * \param [in]        ghost A ghost structure.
 * \return            The total size (in bytes) of the attributes of \a ghost.
 */
size_t              t8_cmesh_trees_ghost_attribute_size (t8_cghost_t ghost);

/* TODO: Currently there is a bug that forces us to give each tree an attribute */
/* TODO: this uses char * and cmesh_set_attribute uses void *. Unify! */
/* attr_tree_index is index of attr in tree's attribute array.
 * We assume that the attributes are already sorted! */
void                t8_cmesh_trees_add_attribute (t8_cmesh_trees_t trees,
                                                  int proc,
                                                  t8_stash_attribute_struct_t
                                                  * attr, t8_locidx_t tree_id,
                                                  size_t index);

/** Return the number of parts of a trees structure.
 * \param [in]        trees The trees structure.
 * \return            The number of parts in \a trees.
 */
size_t              t8_cmesh_trees_get_numproc (t8_cmesh_trees_t trees);

/* TODO: To fit to the interface a trees struct is given as parameter here,
 *       however we could just take the one associated to the cmesh given.*/
/** Print the trees,ghosts and their neighbors in ASCII format t stdout.
 * This function is used for debugging purposes.
 * \param [in]      cmesh A coarse mesh structure that must be committed.
 * \param [in]      trees The trees structure of \a cmesh.
 */
void                t8_cmesh_trees_print (t8_cmesh_t cmesh,
                                          t8_cmesh_trees_t trees);

/** Brodcast an existing valid trees structure from a root rank to
 * all other ranks.
 * The trees structure must belong to cmeshes whose meta_information is
 * already set. \ref t8_cmesh_bcast.
 * \param [in]      cmesh_in    On \a root a committed, replicated cmesh.
 *                              On the other ranks an initialized cmesh with
 *                              the same number of trees as on \a root.
 * \param [in]      root        The rank that broadcasts \a cmesh_in to all
 *                              other ranks.
 * \param [in]      comm        MPI communicator to use.
 */
void                t8_cmesh_trees_bcast (t8_cmesh_t cmesh_in, int root,
                                          sc_MPI_Comm comm);

/** Check whether the face connection of a trees structure are consistent.
 * That is if tree1 lists tree2 as neighbor at face i with ttf entries (or,face j),
 * then tree2 must list tree1 as neighbor at face j with ttf entries (or, face i).
 * \param[in]       cmesh A cmesh structure to be checked.
 * \param[in]       trees The cmesh's trees struct.
 * \return          True if the face connections are consistent,
 *                  False if not.
 */
int                 t8_cmesh_trees_is_face_consistend (t8_cmesh_t cmesh,
                                                       t8_cmesh_trees_t
                                                       trees);

int                 t8_cmesh_trees_is_equal (t8_cmesh_t cmesh,
                                             t8_cmesh_trees_t trees_a,
                                             t8_cmesh_trees_t trees_b);

/** Free all memory allocated with a trees structure.
 *  This means that all coarse trees and ghosts, their face neighbor entries
 *  and attributes and the additional structures of trees are freed.
 * \param [in,out]  trees The tree structure to be destroyed. Set to NULL on output.
 */
void                t8_cmesh_trees_destroy (t8_cmesh_trees_t * trees);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_PART_TREE_H */
