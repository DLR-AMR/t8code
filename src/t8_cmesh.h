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
#include <t8_data/t8_shmem.h>
#include <t8_cmesh/t8_cmesh_save.h>
#include <t8_element.h>

/* TODO: If including eclass were just for the cmesh_new routines, we should
 *       move them into a different file.
 *       However, when specifying the parent-child order in cmesh_reorder,
 *       we might keep the eclass interface for virtual functions.
 *       Actually, we need eclass in the type definition in cmesh.c.
 *       So we might as well use tree-related virtual functions there too.
 */
#include <t8_eclass.h>

/* TODO: See above comment, when moving cmesh_new these get moved too. */
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>

/* TODO: make it legal to call cmesh_set functions multiple times,
 *       just overwrite the previous setting if no inconsistency can occur.
 *       edit: This should be achieved now.
 */

typedef struct t8_cmesh *t8_cmesh_t;
typedef struct t8_ctree *t8_ctree_t;
typedef struct t8_cghost *t8_cghost_t;

T8_EXTERN_C_BEGIN ();

/** Create a new cmesh with reference count one.
 * This cmesh needs to be specialized with the t8_cmesh_set_* calls.
 * Then it needs to be set up with \ref t8_cmesh_commit.
 * \param [in,out] pcmesh       On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new cmesh.
 */
void                t8_cmesh_init (t8_cmesh_t * pcmesh);

/** Check whether a cmesh is not NULL, initialized and not committed.
 * In addition, it asserts that the cmesh is consistent as much as possible.
 * \param [in] cmesh            This cmesh is examined.  May be NULL.
 * \return                      True if cmesh is not NULL,
 *                              \ref t8_cmesh_init has been called on it,
 *                              but not \ref t8_cmesh_commit.
 *                              False otherwise.
 */
int                 t8_cmesh_is_initialized (t8_cmesh_t cmesh);

/** Check whether a cmesh is not NULL, initialized and committed.
 * In addition, it asserts that the cmesh is consistent as much as possible.
 * \param [in] cmesh            This cmesh is examined.  May be NULL.
 * \return                      True if cmesh is not NULL and
 *                              \ref t8_cmesh_init has been called on it
 *                              as well as \ref t8_cmesh_commit.
 *                              False otherwise.
 */
int                 t8_cmesh_is_committed (t8_cmesh_t cmesh);

#ifdef T8_ENABLE_DEBUG
/** After a cmesh is committed, check whether all trees in a cmesh do have positive volume.
 * Returns true if all trees have positive volume.
 * \param [in]  cmesh           This cmesh is examined. May be NULL.
 * \return                      True if \a cmesh is not NULL and all trees for
 *                              which \ref  t8_cmesh_set_tree_vertices
 *                              was called, do have positive geometric volume.
 *                              False otherwise.
 */
int                 t8_cmesh_no_negative_volume (t8_cmesh_t cmesh);
#endif

/** Given a set of vertex coordinates for a tree of a given eclass.
 * Query whether the geometric volume of the tree with this coordinates
 * would be negative.
 * \param [in]  eclass          The eclass of a tree.
 * \param [in]  vertices        The coordinates of the tree's vertices.
 * \param [in]  num_vertices    The number of vertices. \a vertices must hold
 *                              3 * \a num_vertices many doubles.
 *                              \a num_vertices must match \ref t8_eclass_num_vertices[\a eclass]
 * \return                      True if the geometric volume describe by \a vertices is negative.
 *                              Fals otherwise.
 * Returns true if a tree of the given eclass with the given vertex
 * coordinates does have negative volume.
 */
/* TODO: write a test for this function */
int                 t8_cmesh_tree_vertices_negative_volume (t8_eclass_t
                                                            eclass,
                                                            double *vertices,
                                                            int num_vertices);

/* TODO: Currently it is not possible to destroy set_from before
 *       cmesh is destroyed. */
/** This function sets a cmesh to be derived from.
 * The default is to create a cmesh standalone by specifying all data manually.
 * A coarse mesh can also be constructed by deriving it from an existing one.
 * The derivation from another cmesh may optionally be combined with a
 * repartition or uniform refinement of each tree.
 * This function overrides a previously set cmesh to be derived from.
 * \param [in,out] cmesh        Must be initialized, but not committed.
 *                              May even be NULL to revert to standalone.
 * \param [in,out] set_from     Reference counter on this cmesh is bumped.
 *                              It will be unbumped by \ref t8_cmesh_commit,
 *                              after which \a from is no longer remembered.
 *                              Other than that the from object is not changed.
 */
void                t8_cmesh_set_derive (t8_cmesh_t cmesh,
                                         t8_cmesh_t set_from);

/** Allocate a shared memory array to store the tree offsets of a cmesh.
 * \param [in]      mpisize The number of processes.
 * \param [in]      comm    The MPI communicator to use. Its mpisize must match \a mpisize.
 *                  The shared memory type must have been set. Best practice would be
 *                  calling \ref sc_shmem_set_type (comm, T8_SHMEM_BEST_TYPE).
 * \return          A t8_shmem_array struct that stores \a mpisize + 1 t8_gloidx_t entries.
 * \see t8_shmem.h
 */
t8_shmem_array_t    t8_cmesh_alloc_offsets (int mpisize, sc_MPI_Comm comm);

/** Declare if the cmesh is understood as a partitioned cmesh and specify
 * the processor local tree range.
 * This function should be preferred over \ref t8_cmesh_set_partition_offsets
 * when the cmesh is not derived from another cmesh.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
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
 *                             -1: Co not change the face_knowledge level but keep any
 *                                 previously set ones. (Possibly by a previous call to \ref t8_cmesh_set_partition_range)
 * \param [in]     first_local_tree The global index ID of the first tree on this process.
 *                                  If this tree is also the last tree on the previous process,
 *                                  then the argument must be -ID - 1.
 * \param [in]     last_local_tree  The global index of the last tree on this process.
 *                                  If this process should be empty then \a last_local_tree
 *                                  must be strictly smaller than \a first_local_tree.
 *
 * \see t8_cmesh_set_partition_offset \see t8_cmesh_set_partition_uniform
 * \note A value of \a set_face_knowledge other than -1 or 3 is not yet supported.
 */
void                t8_cmesh_set_partition_range (t8_cmesh_t cmesh,
                                                  int set_face_knowledge,
                                                  t8_gloidx_t
                                                  first_local_tree,
                                                  t8_gloidx_t
                                                  last_local_tree);

/* TODO: It is currently not possible to call this function for a non derived
 *       cmesh. Investigate. */
/** Declare if the cmesh is understood as a partitioned cmesh and specify
 * the first local tree for each process.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit.
 * If instead \ref t8_cmesh_set_partition_range was called and the cmesh is
 * derived then the offset array is constructed during commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in] tree_offsets    An array of global tree_id offsets
 *                             for each process can be specified here.
 *                             TODO: document flag for shared trees.
 */
void                t8_cmesh_set_partition_offsets (t8_cmesh_t cmesh,
                                                    t8_shmem_array_t
                                                    tree_offsets);

/** Declare if the cmesh is understood as a partitioned cmesh where the partition
 * table is derived from an assumed uniform refinement of a given level.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     element_level The refinement_level.
 * \param [in]     ts           The element scheme describing the refinement pattern.
 *                              We take ownership. This can be prevented by
 *                              referencing \b ts before calling this function.
 */
void                t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh,
                                                    int element_level,
                                                    t8_scheme_cxx_t * ts);

/* TODO: This function is no longer needed.  Scavenge documentation if helpful. */
#if 0
/* TODO: Currently cmesh_from needs to be partitioned as well.
 *       Change partition function such that it also accepts replicated cmesh_from */
/** Set a cmesh to be partitioned from a second cmesh.
 *  This function can be used instead of \ref t8_cmesh_set_partition.
 *  There a two modes: Either a level is specified, than the new cmesh is partitioned
 *  according to an assumed uniform refinement of the old cmesh,
 *  or an array of tree offsets for each process is specified.
 *  In the latter case each process will get the local trees given by his offsets.
 *  For specification of the offset array see \ref t8_cmesh_types.h.
 * \param [in,out] cmesh       The cmesh to be partitioned.
 * \param [in] cmesh_from      The cmesh to start with.
 * \param [in] level           If >= 0 a uniform refinement of this level is taken
 *                             as reference for the partitioning.
 * \param [in] tree_offsets    If level < 0 then an array of global tree_id offsets
 *                             for each process can be specified here.
 *                             TODO: document flag for shared trees.
 */
void                t8_cmesh_set_partition_given (t8_cmesh_t cmesh,
                                                  t8_gloidx_t * tree_offsets);
#endif

/** Refine the cmesh to a given level.
 * Thus split each tree into x^level subtrees
 * TODO: implement */
/* If level = 0  then no refinement is performed */
void                t8_cmesh_set_refine (t8_cmesh_t cmesh, int level,
                                         t8_scheme_cxx_t * scheme);

/** Set the dimension of a cmesh. If any tree is inserted to the cmesh
 * via \a t8_cmesh_set_tree_class, then the dimension is set automatically
 * to that of the inserted tree.
 * However, if the cmesh is constructed partitioned and the part on this process
 * is empty, it is neccessary to set the dimension by hand.
 * \param [in,out]  cmesh The cmesh to be updated.
 * \param [in]      dim   The dimension to be set. Must satisfy 0 <= dim <= 3.
 * The cmesh must not be committed before calling this function.
 */
void                t8_cmesh_set_dimension (t8_cmesh_t cmesh, int dim);

/** Set the class of a tree in the cmesh.
 * It is not allowed to call this function after \ref t8_cmesh_commit.
 * It is not allowed to call this function multiple times for the same tree.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree_id      The global number of the tree.
 * \param [in]     tree_class   The element class of this tree.
 */
void                t8_cmesh_set_tree_class (t8_cmesh_t cmesh,
                                             t8_gloidx_t gtree_id,
                                             t8_eclass_t tree_class);

/** Store an attribute at a tree in a cmesh.
 *  Attributes can be arbitrary data that is copied to an internal storage
 *  associated to the tree.
 *  Each application can set multiple attributes and attributes are distinguished
 *  by an interger key, where each application can use any integer as key.
 *  TODO: What to do if attribute exists already?
 *        update: Just replace the existing attribute. Our philosophy is that
 *                it is legal to call set functions multiple times.
 *
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
                                            t8_gloidx_t gtree_id,
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
void                t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t gtree1,
                                       t8_gloidx_t gtree2, int face1,
                                       int face2, int orientation);

/** Enable or disable profiling for a cmesh. If profiling is enabled, runtimes
 * and statistics are collected during cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     set_profiling If true, profiling will be enabled, if false
 *                              disabled.
 *
 * Profiling is disabled by default.
 * The cmesh must not be committed before calling this function.
 * \see t8_cmesh_print_profile
 */
void                t8_cmesh_set_profiling (t8_cmesh_t cmesh,
                                            int set_profiling);

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

#ifdef T8_WITH_METIS
/* TODO: document this. */
/* TODO: think about making this a pre-commit set_reorder function. */
void                t8_cmesh_reorder (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/* TODO: think about a sensible interface for a parmetis reordering. */
#endif

/** After allocating and adding properties to a cmesh, finish its construction.
 * TODO: this function is MPI collective.
 * \param [in,out] cmesh        Must be created with \ref t8_cmesh_init
 *                              (TODO: or bcast) and
 *                              specialized with t8_cmesh_set_* calls first (?).
 */
void                t8_cmesh_commit (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/* TODO: Document */
int                 t8_cmesh_save (t8_cmesh_t cmesh, const char *fileprefix);

/* TODO: Document */
t8_cmesh_t          t8_cmesh_load (const char *filename, sc_MPI_Comm comm);

/* TODO: Document */
/* procs_per_node is only relevant in mode==JUQUEEN.
 *  num_files = 1 => replicated cmesh is constructed */
t8_cmesh_t          t8_cmesh_load_and_distribute (const char *fileprefix,
                                                  int num_files,
                                                  sc_MPI_Comm comm,
                                                  t8_load_mode_t mode,
                                                  int procs_per_node);

/** Check whether a given MPI communicator assigns the same rank and mpisize
  * as stored in a cmesh.
  * \param [in] cmesh       The cmesh to be considered.
  * \param [in] comm        A MPI communicator.
  * \return                 True if mpirank and mpisize from \a comm are the same as
  *                         the values stored in \a cmesh.
  *                         False otherwise.
  * \a cmesh must be committed before calling this function.
  * */
int                 t8_cmesh_comm_is_valid (t8_cmesh_t cmesh,
                                            sc_MPI_Comm comm);

/** Query whether a committed cmesh is partitioned or replicated.
 * \param [in] cmesh       A committed cmesh.
 * \return                 True if \a cmesh is partitioned.
 *                         False otherwise.
 * \a cmesh must be committed before calling this function.
 */
int                 t8_cmesh_is_partitioned (t8_cmesh_t cmesh);

/** Return the global number of trees in a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of trees associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_gloidx_t         t8_cmesh_get_num_trees (t8_cmesh_t cmesh);

/** Return the number of local trees of a cmesh.
 *  If the cmesh is not partitioned this is equivalent to \ref t8_cmesh_get_num_trees.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of local trees of the cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t         t8_cmesh_get_num_local_trees (t8_cmesh_t cmesh);

/** Return the number of ghost trees of a cmesh.
 *  If the cmesh is not partitioned this is equivalent to \ref t8_cmesh_get_num_trees.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of ghost trees of the cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t         t8_cmesh_get_num_ghosts (t8_cmesh_t cmesh);

/** Return the global index of the first local tree of a cmesh.
 * If the cmesh is not partitioned this is allways 0.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The global id of the first local tree in cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_gloidx_t         t8_cmesh_get_first_treeid (t8_cmesh_t cmesh);

/** Query whether a given t8_locidx_t belongs to a local tree of a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     An (possible) tree index.
 * \return                 True if \a ltreeid matches the range of local trees of \a cmesh.
 *                         False if not.
 * \a cmesh must be committed before calling this function.
 */
int                 t8_cmesh_treeid_is_local_tree (const t8_cmesh_t cmesh,
                                                   const t8_locidx_t ltreeid);

/** Query whether a given t8_locidx_t belongs to a ghost of a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     An (possible) ghost index.
 * \return                 True if \a ltreeid matches the range of ghost trees of \a cmesh.
 *                         False if not.
 * \a cmesh must be committed before calling this function.
 */
int                 t8_cmesh_treeid_is_ghost (const t8_cmesh_t cmesh,
                                              const t8_locidx_t ltreeid);

/** Given a local tree id that belongs to a ghost, return the index of the ghost.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     The local id of a ghost, satisfying \ref t8_cmesh_treeid_is_ghost,
 *                         thus num_trees <= \a ltreeid < num_trees + num_ghosts
 * \return                 The index of the ghost whithin all ghosts, thus an index
 *                         0 <= index < num_ghosts
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t         t8_cmesh_ltreeid_to_ghostid (const t8_cmesh_t cmesh,
                                                 const t8_locidx_t ltreeid);

/* TODO: Replace this iterator with a new one that does not need the
 *        treeid to be part of the ctree struct */
/* TODO: should this and the next function be part of the interface? */
/** Return a pointer to the first local tree in a cmesh.
 * \param [in]     cmesh        The cmesh to be queried.
 * \return                      A pointer to the first local tree in \a cmesh.
 *                              If \a cmesh has no local trees, NULL is returned.
 * \a cmesh must be committed before calling this function.
 */
t8_ctree_t          t8_cmesh_get_first_tree (t8_cmesh_t cmesh);

/* TODO: should this function behave like first_tree if tree argument is NULL? */
/** Given a local tree in a cmesh return a pointer to the next local tree.
 * \param [in]      cmesh       The cmesh to be queried.
 * \param [in]      tree        A local tree in \a cmesh.
 * \return                      A pointer to the next local tree in \a cmesh
 *                              after \a tree. If no such tree exists, NULL is
 *                              returned.
 * * \a cmesh must be committed before calling this function.
 * TODO: If we run over tree numbers only, don't use ctree_t in API if possible.
 */
t8_ctree_t          t8_cmesh_get_next_tree (t8_cmesh_t cmesh,
                                            t8_ctree_t tree);

/** Return a pointer to a given local tree.
 * \param [in]     cmesh        The cmesh to be queried.
 * \param [in]     ltree_id     The local id of the tree that is asked for.
 * \return                      A pointer to tree in \a cmesh with local
 *                              id \a ltree_id.
 * The cmesh must have at least \a ltree_id + 1 local trees when
 * calling this function.
 * \a cmesh must be committed before calling this function.
 */
t8_ctree_t          t8_cmesh_get_tree (t8_cmesh_t cmesh,
                                       t8_locidx_t ltree_id);

/** Return the eclass of a given local tree.
 * TODO: Should we refer to indices or consequently use ctree_t?
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    tree_id       The local id of the tree whose eclass will be returned.
 * \return                      The eclass of the given tree.
 * TODO: Call tree ids ltree_id or gtree_id etc. instead of tree_id.
 * \a cmesh must be committed before calling this function.
 */
t8_eclass_t         t8_cmesh_get_tree_class (t8_cmesh_t cmesh,
                                             t8_locidx_t ltree_id);

/** Query whether a face of a local tree or ghost is at the domain boundary.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    ltree_id       The local id of a tree.
 * \param [in]    face          The number of a face of the tree.
 * \return                      True if the face is at the domain boundary.
 *                              False otherwise.
 * \a cmesh must be committed before calling this function.
 */
int                 t8_cmesh_tree_face_is_boundary (t8_cmesh_t cmesh,
                                                    t8_locidx_t ltree_id,
                                                    int face);

/** Return the eclass of a given local ghost.
 * TODO: Should we refer to indices or consequently use cghost_t?
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    ghost_id      The local id of the ghost whose eclass will be returned.
 *                              0 <= \a tree_id < cmesh.num_ghosts.
 * \return                      The eclass of the given ghost.
 * \a cmesh must be committed before calling this function.
 */
t8_eclass_t         t8_cmesh_get_ghost_class (t8_cmesh_t cmesh,
                                              t8_locidx_t lghost_id);

/** Return the global id of a given local tree or ghost.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    local_id      The local id of a tree or a ghost.
 *                              If \a local_id < cmesh.num_local_trees then it is
 *                              a tree, otherwise a ghost.
 * \return                      The global id of the tree/ghost.
 */
t8_gloidx_t         t8_cmesh_get_global_id (t8_cmesh_t cmesh,
                                            t8_locidx_t local_id);

/** Return the local id of a give global tree.
 * \param [in]    cmesh         The cmesh.
 * \param [in]    global_id     A global tree id.
 * \return                      Either a value l 0 <= \a l < num_local_trees
 *                              if \a global_id corresponds to a local tree,
 *                              or num_local_trees <= \a l < num_local_trees
 *                              + num_ghosts
 *                              if \a global_id corresponds to a ghost trees,
 *                              or negative if \a global_id neither matches a local
 *                              nor a ghost tree.
 */
t8_locidx_t         t8_cmesh_get_local_id (t8_cmesh_t cmesh,
                                           t8_gloidx_t global_id);

/** Given a local tree id and a face number, get information about the face neighbor tree.
 * \param [in]      cmesh     The cmesh to be considered.
 * \param [in]      ltreeid   The local id of a tree or a ghost.
 * \param [in]      face      A face number of the tree/ghost.
 * \param [out]     dual_face If not NULL, the face number of the neighbor tree at this connection.
 * \param [out]     orientation If not NULL, the face orientation of the connection.
 * \return                    If non-negative: The local id of the neighbor tree or ghost.
 *                            If negative: There is no neighbor across this face. \a dual_face and
 *                            \a orientation remain unchanged.
 * \note If \a ltreeid is a ghost and it has a neighbor which is neither a local tree or ghost,
 *       then the return value will be negative.
 *       Thus, a negative return value does not necessarily mean that this is a domain boundary.
 *       To find out whether a tree is a domain boundary or not \see t8_cmesh_tree_face_is_boundary.
 */
t8_locidx_t         t8_cmesh_get_face_neighbor (const t8_cmesh_t cmesh,
                                                const t8_locidx_t ltreeid,
                                                const int face,
                                                int *dual_face,
                                                int *orientation);

/** Print the collected statistics from a cmesh profile.
 * \param [in]    cmesh         The cmesh.
 *
 * \a cmesh must be committed before calling this function.
 * \see t8_cmesh_set_profiling
 */
void                t8_cmesh_print_profile (t8_cmesh_t cmesh);

/** Return a pointer to the vertex coordinates of a tree.
 * \param [in]    cmesh         The cmesh.
 * \param [in]    ltreeid       The id of a loca tree.
 * \return    If stored, a pointer to the vertex coordinates of \a tree.
 *            If no coordinates for this tree are found, NULL.
 */
double             *t8_cmesh_get_tree_vertices (t8_cmesh_t cmesh,
                                                t8_locidx_t ltreeid);

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
                                            t8_locidx_t ltree_id);

/** Return the shared memory array storing the partition table of
 * a partitioned cmesh.
 * \param [in]      cmesh       The cmesh.
 * \return                      The partition array.
 *                              NULL if the cmesh is not partitioned or
 *                              the partition array is not stored in \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_shmem_array_t    t8_cmesh_get_partition_table (t8_cmesh_t cmesh);

/* TODO: remove get_ when there is no risk of confusion? Convention?
 *       Update: use get throughout for access functions that do not change the object.
 * */

/** Calculate the section of a uniform forest for the current rank.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    level         The uniform refinement level to be created.
 * \param [in]    ts            The element scheme for which to compute the bounds.
 * \param [out]   first_local_tree  The first tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_begin The global index of the first element belonging to the calling processor. Not computed if NULL.
 * \param [out]   last_local_tree  The last tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_end The global index of the first element that does not belonging to
 *                                  the calling processor anymore. Not computed if NULL.
 * \param [out]   first_tree_shared If not NULL, 1 or 0 is stored here depending on whether \a first_local_tree is the
 *                                 same as \a last_local_tree on the next process.
 * \a cmesh must be committed before calling this function. *
 */
void                t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                                             t8_scheme_cxx_t * ts,
                                             t8_gloidx_t * first_local_tree,
                                             t8_gloidx_t *
                                             child_in_tree_begin,
                                             t8_gloidx_t * last_local_tree,
                                             t8_gloidx_t * child_in_tree_end,
                                             int8_t * first_tree_shared);

/** Increase the reference counter of a cmesh.
 * \param [in,out] cmesh        On input, this cmesh must exist with positive
 *                              reference count.  It may be in any state.
 */
void                t8_cmesh_ref (t8_cmesh_t cmesh);

/** Decrease the reference counter of a cmesh.
 * If the counter reaches zero, this cmesh is destroyed.
 * See also \ref t8_cmesh_destroy, which is to be preferred when it is
 * known that the last reference to a cmesh is deleted.
 * \param [in,out] pcmesh       On input, the cmesh pointed to must exist
 *                              with positive reference count.  It may be in
 *                              any state.  If the reference count reaches
 *                              zero, the cmesh is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the cmesh is not modified in other ways.
 */
void                t8_cmesh_unref (t8_cmesh_t * pcmesh);

/** Verify that a coarse mesh has only one reference left and destroy it.
 * This function is preferred over \ref t8_cmesh_unref when it is known
 * that the last reference is to be deleted.
 * \param [in,out]  pcmesh      This cmesh must have a reference count of one.
 *                              It can be in any state (committed or not).
 *                              Then it effectively calls \ref t8_cmesh_unref.
 * \param [in]      comm        A mpi communicator that is valid with \a cmesh.
 */
void                t8_cmesh_destroy (t8_cmesh_t * pcmesh);

/* Functions for construcing complete and committed cmeshes */

/** Constructs a cmesh from a given p4est_connectivity structure.
 * \param[in]       conn       The p4est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p4est (p4est_connectivity_t * conn,
                                             sc_MPI_Comm comm,
                                             int do_partition);

/** Constructs a cmesh from a given p8est_connectivity structure.
 * \param[in]       conn       The p8est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_dup     Flag whether the communicator shall be duplicated or not.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 */
t8_cmesh_t          t8_cmesh_new_from_p8est (p8est_connectivity_t * conn,
                                             sc_MPI_Comm comm,
                                             int do_partition);

/* TODO: it could possibly be a problem that we do not set the dimension of
 * the cmesh. This could i.e. be difficult when we combine an empty cmesh with
 * a non-empty one. */
/** Construct a cmesh that has no trees. We do not know a special use case,
 * this function is merely for debugging and to show the possibility.
 * \param [in]      comm       mpi communicator to be used with the new cmesh.
 * \param [in]      do_partition Flag whether the cmesh should be partitioned or not.
 * \return                     A committed t8_cmesh structure that has no trees.
 */
t8_cmesh_t          t8_cmesh_new_empty (sc_MPI_Comm comm, int do_partition);

/** Constructs a cmesh that consists only of one tree of a given element class.
 * \param [in]      eclass     The element class.
 * \param [in]      comm       mpi communicator to be used with the new cmesh.
 * \param [in]      do_dup     Flag whether the communicator shall be duplicated or not.
 * \return          A committed t8_cmesh structure with one tree of class \a eclass.
 */
t8_cmesh_t          t8_cmesh_new_from_class (t8_eclass_t eclass,
                                             sc_MPI_Comm comm);

t8_cmesh_t          t8_cmesh_new_testhybrid (sc_MPI_Comm comm);

/** Construct a hypercube forest from one primitive tree class.
 * \param [in] eclass       This element class determines the dimension and
 *                          the number of trees needed to construct a cube.
 * \param [in] comm         The mpi communicator to be used.
 * \param [in] do_bcast     If this flag is nonzero the cmesh is only constructed
 *                          on processor 0 and then broadcasted to the other
 *                          processors in \a comm.
 *                          TODO: this parameter will be moved to internal.
 * \param [in] do_partition Create a partitioned cmesh.
 * \param [in] periodic     If true, the coarse mesh will be periodic in each direction.
 *                          Not possible with \a eclass pyramid.
 * TODO: Add periodic flags for each dimension.
 */
t8_cmesh_t          t8_cmesh_new_hypercube (t8_eclass_t eclass,
                                            sc_MPI_Comm comm,
                                            int do_bcast, int do_partition,
                                            int periodic);

/** Hybercube with 6 Tets, 6 Prism, 4 Hex. */
/* TODO: Document */
t8_cmesh_t          t8_cmesh_new_hypercube_hybrid (int dim, sc_MPI_Comm comm,
                                                   int do_partition,
                                                   int periodic);

/** Construct a unit interval/square/cube coarse mesh that is periodic in each direction.
 * Element class?
 * Hypercube?
 * TODO: redundant, remove.
 * \param [in] comm         The mpi communicator to use.
 * \param [in] dim          The dimension of the forest, 1, 2 or 3.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic (sc_MPI_Comm comm, int dim);

/** Construct a unit square of two triangles that is periodic in x and y.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_tri (sc_MPI_Comm comm);

/** Construct a unit square of two quads and four triangles that is periodic in x and y.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_hybrid (sc_MPI_Comm comm);

/** Construct a unit interval coarse mesh that consists of 3 trees and is
 * periodic.
 * \param [in] comm         The mpi communicator to use.
 * \return                  A valid cmesh, as is _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_periodic_line_more_trees (sc_MPI_Comm comm);

/** Construct a mesh consisting of a given number of same type trees.
 * \param [in] eclass       This element class determines the dimension and
 *                          the type trees used.
 * \param [in] num_trees    The number of trees to use.
 * \param [in] comm         The MPI_Communicator used to commit the cmesh.
 * \return                  A valid cmesh, as if _init and _commit had been called.
 */
t8_cmesh_t          t8_cmesh_new_bigmesh (t8_eclass_t eclass, int num_trees,
                                          sc_MPI_Comm comm);

/** Construct a forest of three connected askew lines
  * \param [in] comm         The mpi communicator to use.
  * \return                  A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_line_zigzag (sc_MPI_Comm comm);

/** Construct a forest of num_of_prisms connected prism, all with one edge in 0,
  * except for num_of_prisms = 2, then the return is the hypercube mesh
  * \param [in] comm        The mpi communicator to use.
  * \param [in] num_of_prisms The number of prisms to be used.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_cake (sc_MPI_Comm comm,
                                             int num_of_prisms);

/** Construct a single deformed prism
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh; as if _init and _commit had been called.*/
t8_cmesh_t          t8_cmesh_new_prism_deformed (sc_MPI_Comm comm);

/** Construct a forest of six connected noncannoical oriented prisms
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_cake_funny_oriented (sc_MPI_Comm comm);

/** Construct a forest of six connected noncannoical oriented prisms
  * \param [in] comm        The mpi communicator to use.
  * \return                 A valid cmesh, as if _init and _commit had been called.
  */
t8_cmesh_t          t8_cmesh_new_prism_geometry (sc_MPI_Comm comm);

/** Create a partitoned cmesh of quads whose local trees are given by an
 * num_x by num_y brick connectivity from p4est
 * or a num_x by num_y by num_z brick connectivity from p8est.
 * num_x and num_y and num_z can be different for different MPI ranks.
 * \param [in] num_x       The number of trees in x direction for this rank. Must be >= 0.
 * \param [in] num_y       The number of trees in y direction for this rank. Must be >= 0.
 * \param [in] num_y       The number of trees in z direction for this rank. Must be >= 0.
 *                         If nonzero, the cmesh is 3 dimensional.
 * \param [in] x_periodic  If nonzero, the local brick connectivity is periodic in x direction.
 * \param [in] y_periodic  If nonzero, the local brick connectivity is periodic in y direction.
 * \param [in] y_periodic  If nonzero and \a num_z > 0, the local brick connectivity is periodic in z direction.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and partitioned cmesh. The process local trees
 *                         form a \a num_x by \a num_y (by \a num_z) brick.
 * It is possible for num_x or num_y to be set to zero. In this case the local part
 * of the cmesh will be empty.
 * If num_z is set to zero, the cmesh is 2 dimensional.
 */
t8_cmesh_t          t8_cmesh_new_disjoint_bricks (t8_gloidx_t num_x,
                                                  t8_gloidx_t num_y,
                                                  t8_gloidx_t num_z,
                                                  int x_periodic,
                                                  int y_periodic,
                                                  int z_periodic,
                                                  sc_MPI_Comm comm);

/** Construct a tetrahedral cmesh that has all possible face to face
 * connections and orientations.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated cmesh of 24 tetrahedron trees
 *                         in which each (face -> face, orientation) face connection
 *                         is set at least once.
 *                         Note that most faces in this cmesh are boundary faces.
 */
t8_cmesh_t          t8_cmesh_new_tet_orientation_test (sc_MPI_Comm comm);

/** Construct a hybrid cmesh with 2 tets, 2 prism, 1 hex.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated hybrid cmesh of 5 trees.
 */
t8_cmesh_t          t8_cmesh_new_hybrid_gate (sc_MPI_Comm comm);

/** Construct a hybrid cmesh with 2 tets, 2 prism, 1 hex and all are deformed.
 * This cmesh is used for testing and debugging.
 * \param [in] comm        The MPI communicator used to commit the cmesh.
 * \return                 A committed and replicated hybrid cmesh of 5 trees.
 */
t8_cmesh_t          t8_cmesh_new_hybrid_gate_deformed (sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_H */
