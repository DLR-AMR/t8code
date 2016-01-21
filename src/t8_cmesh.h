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

/** \file t8_cmesh.h
 * We define the coarse mesh of trees in this file.
 */

#ifndef T8_CMESH_H
#define T8_CMESH_H

#include <t8.h>
#include <t8_eclass.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>

/* TODO: do set_mpicomm and figure out dup logic */

typedef struct t8_cmesh *t8_cmesh_t;
typedef struct t8_ctree *t8_ctree_t;
typedef struct t8_cghost *t8_cghost_t;

T8_EXTERN_C_BEGIN ();

/** Create a new cmesh with reference count one.
 * This cmesh needs to be specialized with the t8_cmesh_set_* calls.
 * Then it needs to be set up with \see t8_cmesh_commit.
 * \param [in,out] pcmesh       On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new cmesh.
 */
void                t8_cmesh_init (t8_cmesh_t * pcmesh);

/** Set MPI communicator to use in commiting a new cmesh.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \see t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh whose communicator will be set.
 * \param [in] mpicomm          This MPI communicator must be valid.
 * \param [in] do_dup           If true, the communicator will be duped in
 *                              this creation and whenever another cmesh is
 *                              derived from it in the future.
 *                              If false, no duping takes place at all.
 */
void                t8_cmesh_set_mpicomm (t8_cmesh_t cmesh,
                                          sc_MPI_Comm mpicomm, int do_dup);

/* TODO: think about how set_num_trees and this function play together */
/** Declare if the cmesh is understood as a partitioned cmesh or a
 * replicated cmesh. Replicated (each processor owns the whole mesh) is
 * the default and in this case \ref t8_cmesh_set_partitioned is the same as
 * \ref t8_cmesh_set_num_trees and the values \a first_local_tree and
 * \a set_face_knowledge are ignored.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \see t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     set_partitioned A nonzero value specifies that \a cmesh
 *                              is interpreted as a partitioned mesh.
 * \parma [in]     set_face_knowledge   Several values are possible that define
 *                              how much information is required on face connections,
 *                              specified by \ref t8_cmesh_set_join.
 *                              0: Expect face connection of local trees.
                                1: In addition, expect face connection from
 *                                 ghost trees to local trees.
 *                              2: In addition, expect face connection between
 *                                 ghost trees.
 *                              3: Expect face connection of local and ghost trees.
 *                              Consistency of this requirement is checked on
 *                              \ref t8_cmesh_commit.
 * \param [in]     first_local_tree The global index of the first tree on this process.
 *                                  If \a set_partitioned is zero, must be 0.
 * \param [in]     last_local_tree  The global index of the last tree on this process.
 *                                  If \a set_partitioned is zero, must be
 *                                  \a num_global_trees - 1.
 */
void                t8_cmesh_set_partitioned (t8_cmesh_t cmesh,
                                              int set_partitioned,
                                              int set_face_knowledge,
                                              t8_gloidx_t first_local_tree,
                                              t8_gloidx_t last_local_tree);

/* TODO: document.
 *       if level >= 0 then ignore trees_per_proc
 */
void                t8_cmesh_set_partition_from (t8_cmesh_t cmesh,
                                                 const t8_cmesh_t cmesh_from,
                                                 int level,
                                                 t8_locidx_t * trees_per_proc);

/* TODO: This is actually not part of the interface?
 *       At least it is only used if the cmesh is partitioned.
 */
/** Set the total number of trees for a coarse mesh.
 * It is not allowed to call this function after \ref t8_cmesh_commit.
 * TODO: Clarify that this holds for all _set_ functions.
 * TODO: Explain what _commit does in general as a convention.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     num_trees    The total number of trees in this coarse mesh.
 */
void                t8_cmesh_set_num_trees (t8_cmesh_t cmesh,
                                            t8_gloidx_t num_trees);

/** Set the class of a tree in the cmesh.
 * It is not allowed to call this function after \ref t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree_id      The global number of the tree.
 * \param [in]     tree_class   The element class of this tree.
 */
void                t8_cmesh_set_tree_class (t8_cmesh_t cmesh,
                                             t8_gloidx_t tree_id,
                                             t8_eclass_t tree_class);

/** Store an attribute at a tree in a cmesh.
 *  Attributes can be arbitrary data that is copied to an internal storage
 *  associated to the tree.
 *  Each application can set multiple attributes and attributes are distinguished
 *  by an interger key, where each application can use any integer as key.
 *  TODO: What to do if attribute exists already?
 * \param [in, out] cmesh       The cmesh to be updated.
 * \param [in]      tree_id     The global id of the tree.
 * \param [in]      package_id  Unique identifier of a valid software package. \see sc_package_register
 * \param [in]      key         An integer key used to identify this attribute under all
 *                              attributes with the same package_id.
 *                              \a key must be a unique value for this tree and package_id.
 * \param [in]      data        A pointer to the attribute data.
 * \param [in]      data_size   The number of bytes of the attribute.
 * \param [in]      data_persists This flag can be used to optimize memory. If true
 *                              then t8code assumes that the attribute data is present at the
 *                              memory that \a data points to when \ref t8_cmesh_commit is called
 *                              (This is more memory efficient).
 *                              If the flag is false an internal copy of the data is created
 *                              immediately and this copy is used at commit.
 *                              In both cases a copy of the data is used by t8_code after t8_cmesh_commit.
 */
void                t8_cmesh_set_attribute (t8_cmesh_t cmesh,
                                            t8_gloidx_t tree_id,
                                            int package_id, int key,
                                            void *data, size_t data_size,
                                            int data_persists);

/** Insert a face-connection between two trees in a cmesh.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree1        The tree id of the first of the two trees.
 * \param [in]     tree2        The tree id of the second of the two trees.
 * \param [in]     face1        The face number of the first tree.
 * \param [in]     face2        The face number of the second tree.
 * \param [in]     orientation  Specify how face1 and face2 are oriented to each other
 *                              TODO: orientation needs to be carefully defined
 *                              for all element classes.
 * TODO: document orientation
 */
void                t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t tree1,
                                       t8_gloidx_t tree2, int face1,
                                       int face2, int orientation);

/* returns true if cmesh_a equals cmesh_b */
/* TODO: document
 * collective or serial */
/** Check whether two given cmeshes carry the same information.
 * \param [in]    cmesh_a       The first of the two cmeshes to be checked.
 * \param [in]    cmesh_b       The second of the two cmeshes to be checked.
 * \return                      True if both cmeshes carry the same information,
 *                              false otherwise.
 *                              TODO: define carefully.
 *                              Orders, sequences, equivalences?
 * Currently the attributes of the trees are not compared.
 * This function works on committed and uncommitted cmeshes.
 */
int                 t8_cmesh_is_equal (t8_cmesh_t cmesh_a,
                                       t8_cmesh_t cmesh_b);

/** Broadcast a cmesh structure that exists only on one process to all
 *  processes in the cmesh's communicator.
 *  TODO: Input structure must be replicated, not parallelized.
 *  TODO: Recommend to call this just before commit.  Earlier is thinkable too.
 *  On the other processors, it will be allocated.
 *  It is not allowed to call this function after \ref t8_cmesh_commit.
 *  \param [in] cmesh_in For the root process the cmesh to be broadcast,
 *                      for the other processes it must be NULL.
 *  \param [in] root    The rank of the process that provides the cmesh.
 *  \param [in] comm    The mpi communicator. Must match cmesh's communicator
 *                      on the root process.
 *  \return             For the root process this is a pointer to \a cmesh_in.
 *                      Else, a pointer to a newly allocated cmesh
 *                      structure with the same values as \a conn_in on the
 *                      root process.
 */
t8_cmesh_t          t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root,
                                    sc_MPI_Comm comm);

/** After allocating and adding properties to a cmesh, finish its construction.
 * TODO: this function is MPI collective.
 * \param [in,out] cmesh        Must be created with \ref t8_cmesh_init
 *                              (TODO: or bcast) and
 *                              specialized with t8_cmesh_set_* calls first (?).
 */
void                t8_cmesh_commit (t8_cmesh_t cmesh);

/** Return the MPI communicator of a cmesh.
 * TODO: Move all _get_ functions after _commit.
 * \param [in] cmesh       The cmesh whose communicator will be returned.
 * \param [out] do_dup     This variable is filled with the do_dup entry of \a cmesh.
 * \return                 The MPI communicator associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
sc_MPI_Comm         t8_cmesh_get_mpicomm (t8_cmesh_t cmesh, int *do_dup);

/** Return the global number of trees in a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of trees associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_gloidx_t         t8_cmesh_get_num_trees (t8_cmesh_t cmesh);

/** Return the processor local number of trees in a cmesh.
 * If the cmesh is not partitioned this is the same as
 * the global number of trees.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of trees associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t         t8_cmesh_get_local_num_trees (t8_cmesh_t cmesh);

/* TODO: should this and the next function be part of the interface? */
/** Return a pointer to the first local tree in a cmesh.
 * \param [in]     cmesh        The cmesh to be queried.
 * \return                      A pointer to the first local tree in \a cmesh.
 *                              If \a cmesh has no local trees, NULL is returned.
 * \a cmesh must be committed before calling this function.
 */
t8_ctree_t          t8_cmesh_first_tree (t8_cmesh_t cmesh);

/* TODO: should this function behave like first_tree if tree argument is NULL? */
/** Given a local tree in a cmesh return a pointer to the next local tree.
 * \param [in]      cmesh       The cmesh to be queried.
 * \param [in]      tree        A local tree in \a cmesh.
 * \return                      A pointer to the next local tree in \a cmesh
 *                              after \a tree. If no such tree exists, NULL is
 *                              returned.
 * * \a cmesh must be committed before calling this function.
 */
t8_ctree_t          t8_cmesh_next_tree (t8_cmesh_t cmesh, t8_ctree_t tree);

/** Return the eclass of a given local tree.
 * TODO: Should we refer to indices or consequently use ctree_t?
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    tree_id       The local id of the tree whose eclass will be returned.
 * \return                      The eclass of the given tree.
 * \a cmesh must be committed before calling this function.
 */
t8_eclass_t         t8_cmesh_get_tree_class (t8_cmesh_t cmesh,
                                             t8_locidx_t tree_id);

/** Return the attribute pointer of a tree.
 * \param [in]     cmesh        The cmesh.
 * \param [in]     package_id   The identifier of a valid software package. \see sc_package_register
 * \param [in]     key          A key used to identify the attribute under all
 *                              attributes of this tree with the same \a package_id.
 * \param [in]     tree_id      The local number of the tree.
 * \param [out]    data_size    The size of the attribute in bytes.
 * \return         The attribute pointer of the tree \a tree_id.
 * \a cmesh must be committed before calling this function.
 * \see t8_cmesh_set_attribute
 */
void               *t8_cmesh_get_attribute (t8_cmesh_t cmesh,
                                            int package_id, int key,
                                            t8_locidx_t tree_id,
                                            size_t * data_size);

/* TODO: remove get_ when there is no risk of confusion? Convention? */

/** Calculate the section of a uniform forest for the current rank.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    level         The uniform refinement level to be created.
 * \param [out]   first_local_tree  The first tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_begin The global index of the first element belonging to the calling processor.
 * \param [out]   last_local_tree  The last tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_end The global index of the last element belonging to the calling processor.
 * \a cmesh must be committed before calling this function.
 */
void                t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                                             t8_gloidx_t * first_local_tree,
                                             t8_gloidx_t *
                                             child_in_tree_begin,
                                             t8_gloidx_t * last_local_tree,
                                             t8_gloidx_t * child_in_tree_end);

/** Increase the reference counter of a cmesh.
 * \param [in,out] cmesh        On input, this cmesh must exist with positive
 *                              reference count.  It may be in any state.
 */
void                t8_cmesh_ref (t8_cmesh_t cmesh);

/** Decrease the reference counter of a cmesh.
 * If the counter reaches zero, this cmesh is destroyed.
 * \param [in,out] pcmesh       On input, the cmesh pointed to must exist
 *                              with positive reference count.  It may be in
 *                              any state.  If the reference count reaches
 *                              zero, the cmesh is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the cmesh is not modified in other ways.
 */
void                t8_cmesh_unref (t8_cmesh_t * pcmesh);

/* Functions for construcing complete and committed cmeshes */

/** Constructs a cmesh from a given p4est_connectivity structure.
 *  The constructed cmesh will be replicated.
 * \param[in]       conn       The p4est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_dup     Flag whether the communicator shall be duplicated or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p4est (p4est_connectivity_t * conn,
                                             sc_MPI_Comm comm, int do_dup);

/** Constructs a cmesh from a given p8est_connectivity structure.
 *  The constructed cmesh will be replicated.
 * \param[in]       conn       The p8est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_dup     Flag whether the communicator shall be duplicated or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p8est (p8est_connectivity_t * conn,
                                             sc_MPI_Comm comm, int do_dup);

/** Constructs a cmesh that consists only of one tree of a given element class.
 * \param [in]      eclass     The element class.
 * \param [in]      comm       mpi communicator to be used with the new cmesh.
 * \param [in]      do_dup     Flag whether the communicator shall be duplicated or not.
 * \return          A committed t8_cmesh structure with one tree of class \a eclass.
 */
t8_cmesh_t          t8_cmesh_new_from_class (t8_eclass_t eclass,
                                             sc_MPI_Comm comm, int do_dup);

/** Construct a hypercube forest from one primitive tree class.
 * \param [in] eclass       This element class determines the dimension and
 *                          the number of trees needed to construct a cube.
 * \param [in] comm         The mpi communicator to be used.
 * \param [in] do_dup       Whether \a comm is to be duplicated.
 * \param [in] do_bcast     If this flag is nonzero the cmesh is only constructed
 *                          on processor 0 and then broadcasted to the other
 *                          processors in \a comm.
 *                          TODO: this parameter will be moved to internal.
 * \param [in] do_partition Create a partitioned cmesh.
 * TODO: Add periodic flags for each dimension.
 */
t8_cmesh_t          t8_cmesh_new_hypercube (t8_eclass_t eclass,
                                            sc_MPI_Comm comm, int do_dup,
                                            int do_bcast, int do_partition);

/** Construct a unit interval/square/cube forest that is periodic in each direction.
 * Element class?
 * Hypercube?
 * TODO: redundant, remove.
 * \param [in] comm         The mpi communicator to use.
 * \param [in] do_dup       Whether the mpi communicator is to be duplicated.
 * \param [in] dim          The dimension of the forest, 1, 2 or 3.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic (sc_MPI_Comm comm, int do_dup,
                                           int dim);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_H */
