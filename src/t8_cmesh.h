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
#include <t8_schemes/t8_scheme.h>

/* Forward pointer reference to hidden cmesh implementation.
 * This reference needs to be known by t8_geometry, hence we 
 * put it before the include. */
typedef struct t8_cmesh *t8_cmesh_t;

#include <t8_geometry/t8_geometry.h>

/* TODO: If including eclass were just for the cmesh_new routines, we should
 *       move them into a different file.
 *       However, when specifying the parent-child order in cmesh_reorder,
 *       we might keep the eclass interface for virtual functions.
 *       Actually, we need eclass in the type definition in cmesh.c.
 *       So we might as well use tree-related virtual functions there too.
 */
#include <t8_eclass.h>

/* TODO: make it legal to call cmesh_set functions multiple times,
 *       just overwrite the previous setting if no inconsistency can occur.
 *       edit: This should be achieved now.
 */

/* Forward pointer references to hidden implementations of
 * tree and ghost tree. */
typedef struct t8_ctree *t8_ctree_t;
typedef struct t8_cghost *t8_cghost_t;

T8_EXTERN_C_BEGIN ();

/** Create a new cmesh with reference count one.
 * This cmesh needs to be specialized with the t8_cmesh_set_* calls.
 * Then it needs to be set up with \ref t8_cmesh_commit.
 * \param [in,out] pcmesh       On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new cmesh.
 */
void
t8_cmesh_init (t8_cmesh_t *pcmesh);

/** Allocate a new un-committed cmesh.
 * \return                     A pointer to an un-committed t8_cmesh structure.
 */
t8_cmesh_t
t8_cmesh_new ();

/** Check whether a cmesh is not NULL, initialized and not committed.
 * In addition, it asserts that the cmesh is consistent as much as possible.
 * \param [in] cmesh            This cmesh is examined.  May be NULL.
 * \return                      True if cmesh is not NULL,
 *                              \ref t8_cmesh_init has been called on it,
 *                              but not \ref t8_cmesh_commit.
 *                              False otherwise.
 */
int
t8_cmesh_is_initialized (t8_cmesh_t cmesh);

/** Check whether a cmesh is not NULL, initialized and committed.
 * In addition, it asserts that the cmesh is consistent as much as possible.
 * \param [in] cmesh            This cmesh is examined.  May be NULL.
 * \return                      True if cmesh is not NULL and
 *                              \ref t8_cmesh_init has been called on it
 *                              as well as \ref t8_cmesh_commit.
 *                              False otherwise.
 */
int
t8_cmesh_is_committed (const t8_cmesh_t cmesh);

#if T8_ENABLE_DEBUG
/** Check the geometry of the mesh for validity.
 * \param [in] cmesh            This cmesh is examined.
 * \return                      True if the geometry of the cmesh is valid.
 */
int
t8_cmesh_validate_geometry (const t8_cmesh_t cmesh);
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
 *                              False otherwise.
 * Returns true if a tree of the given eclass with the given vertex
 * coordinates does have negative volume.
 */
/* TODO: write a test for this function */
int
t8_cmesh_tree_vertices_negative_volume (const t8_eclass_t eclass, const double *vertices, const int num_vertices);

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
void
t8_cmesh_set_derive (t8_cmesh_t cmesh, t8_cmesh_t set_from);

/** Allocate a shared memory array to store the tree offsets of a cmesh.
 * \param [in]      mpisize The number of processes.
 * \param [in]      comm    The MPI communicator to use. Its mpisize must match \a mpisize.
 *                  The shared memory type must have been set. Best practice would be
 *                  calling \ref sc_shmem_set_type (comm, T8_SHMEM_BEST_TYPE).
 * \return          A t8_shmem_array struct that stores \a mpisize + 1 t8_gloidx_t entries.
 * \see t8_shmem.h
 */
t8_shmem_array_t
t8_cmesh_alloc_offsets (int mpisize, sc_MPI_Comm comm);

/** Declare if the cmesh is understood as a partitioned cmesh and specify
 * the processor local tree range.
 * This function should be preferred over \ref t8_cmesh_set_partition_offsets
 * when the cmesh is not derived from another cmesh.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     set_face_knowledge   Several values are possible that define
 *                              how much information is required on face connections,
 *                              specified by \ref t8_cmesh_set_join.
 *                              0: Expect face connection of local trees.
 *                              1: In addition, expect face connection from
 *                                 ghost trees to local trees.
 *                              2: In addition, expect face connection between
 *                                 ghost trees.
 *                              3: Expect face connection of local and ghost trees.
 *                              Consistency of this requirement is checked on
 *                              \ref t8_cmesh_commit.
 *                             -1: Do not change the face_knowledge level but keep any
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
void
t8_cmesh_set_partition_range (t8_cmesh_t cmesh, int set_face_knowledge, t8_gloidx_t first_local_tree,
                              t8_gloidx_t last_local_tree);

/** Declare if the cmesh is understood as a partitioned cmesh and specify
 * the first local tree for each process.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit.
 * If instead \ref t8_cmesh_set_partition_range was called and the cmesh is
 * derived then the offset array is constructed during commit.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in] tree_offsets     An array of global tree_id offsets
 *                              for each process can be specified here.
 *                             TODO: document flag for shared trees.
 */
void
t8_cmesh_set_partition_offsets (t8_cmesh_t cmesh, t8_shmem_array_t tree_offsets);

/** Declare if a derived cmesh should be partitioned according to a
 * uniform refinement of a given level for the provided scheme.
 * This call is only valid when the cmesh is not yet committed via a call
 * to \ref t8_cmesh_commit and when the cmesh will be derived.
 * \param [in,out] cmesh          The cmesh to be updated.
 * \param [in]     element_level  The refinement_level.
 * \param [in]     scheme             The element scheme describing the refinement pattern.
 *                                We take ownership. This can be prevented by
 *                                referencing \b scheme before calling this function.
 */
void
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, const int element_level, const t8_scheme_c *scheme);

/** Refine the cmesh to a given level.
 * Thus split each tree into x^level subtrees
 * TODO: implement */
/* If level = 0  then no refinement is performed */
void
t8_cmesh_set_refine (t8_cmesh_t cmesh, const int level, const t8_scheme_c *scheme);

/** Set the dimension of a cmesh. If any tree is inserted to the cmesh
 * via \a t8_cmesh_set_tree_class, then the dimension is set automatically
 * to that of the inserted tree.
 * However, if the cmesh is constructed partitioned and the part on this process
 * is empty, it is necessary to set the dimension by hand.
 * \param [in,out]  cmesh The cmesh to be updated.
 * \param [in]      dim   The dimension to be set. Must satisfy 0 <= dim <= 3.
 * The cmesh must not be committed before calling this function.
 */
void
t8_cmesh_set_dimension (t8_cmesh_t cmesh, int dim);

/** Set the class of a tree in the cmesh.
 * It is not allowed to call this function after \ref t8_cmesh_commit.
 * It is not allowed to call this function multiple times for the same tree.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree_id      The global number of the tree.
 * \param [in]     tree_class   The element class of this tree.
 */
void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, t8_eclass_t tree_class);

/** Store an attribute at a tree in a cmesh.
 *  Attributes can be arbitrary data that is copied to an internal storage
 *  associated to the tree.
 *  Each application can set multiple attributes and attributes are distinguished
 *  by an integer key, where each application can use any integer as key.
 *
 * \param [in, out] cmesh       The cmesh to be updated.
 * \param [in]      gtree_id    The global id of the tree.
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
 * \note If an attribute with the given package_id and key already exists, then it will get overwritten.
 */
void
t8_cmesh_set_attribute (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key, void *data, size_t data_size,
                        int data_persists);

/** Store a string as an attribute at a tree in a cmesh.
 * \param [in, out] cmesh       The cmesh to be updated.
 * \param [in]      gtree_id    The global id of the tree.
 * \param [in]      package_id  Unique identifier of a valid software package. \see sc_package_register
 * \param [in]      key         An integer key used to identify this attribute under all
 *                              attributes with the same package_id.
 *                              \a key must be a unique value for this tree and package_id.
 * \param [in]      string      The string to store as attribute.
 * \note You can also use \ref t8_cmesh_set_attribute, but we recommend using this
 *       specialized function for strings.
 * \note If an attribute with the given package_id and key already exists, then it will get overwritten.
 */
void
t8_cmesh_set_attribute_string (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key, const char *string);

/** Store an array of t8_gloidx_t as an attribute at a tree in a cmesh.
 * \param [in, out] cmesh       The cmesh to be updated.
 * \param [in]      gtree_id    The global id of the tree.
 * \param [in]      package_id  Unique identifier of a valid software package. \see sc_package_register
 * \param [in]      key         An integer key used to identify this attribute under all
 *                              attributes with the same package_id.
 *                              \a key must be a unique value for this tree and package_id.
 * \param [in]      data        The array to store as attribute.
 * \param [in]      data_count  The number of entries in \a data.
 * \param [in]      data_persists This flag can be used to optimize memory. If true
 *                              then t8code assumes that the attribute data is present at the
 *                              memory that \a data points to when \ref t8_cmesh_commit is called
 *                              (This is more memory efficient).
 *                              If the flag is false an internal copy of the data is created
 *                              immediately and this copy is used at commit.
 *                              In both cases a copy of the data is used by t8_code after t8_cmesh_commit.
 * \note You can also use \ref t8_cmesh_set_attribute, but we recommend using this
 *       specialized function for arrays.
 * \note If an attribute with the given package_id and key already exists, then it will get overwritten.
 * \note We do not store the number of data entries \a data_count of the attribute array.
 *       You can keep track of the data count yourself by using another attribute.
 */
void
t8_cmesh_set_attribute_gloidx_array (t8_cmesh_t cmesh, t8_gloidx_t gtree_id, int package_id, int key,
                                     const t8_gloidx_t *data, const size_t data_count, int data_persists);

/** Insert a face-connection between two trees in a cmesh.
 * \param [in,out] cmesh        The cmesh to be updated.
 * \param [in]     tree1        The tree id of the first of the two trees.
 * \param [in]     tree2        The tree id of the second of the two trees.
 * \param [in]     face1        The face number of the first tree.
 * \param [in]     face2        The face number of the second tree.
 * \param [in]     orientation  Specify how face1 and face2 are oriented to each other
 * 
 * \note The orientation is defined as:
 * Let my_face and other_face be the two face numbers of the connecting trees.
 * We chose a main_face from them as follows: Either both trees have the same
 * element class, then the face with the lower face number is the main_face or
 * the trees belong to different classes in which case the face belonging to the
 * tree with the lower class according to the ordering
 * triangle < quad, hex < tet < prism < pyramid, is the main_face.
 * Then face corner 0 of the main_face connects to a face
 * corner k in the other face.  The face orientation is defined as the number k.
 * If the classes are equal and my_face == other_face, treating
 * either of both faces as the main_face leads to the same result.
 * See https://arxiv.org/pdf/1611.02929.pdf for more details.
 */
void
t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t gtree1, t8_gloidx_t gtree2, int face1, int face2, int orientation);

/** Enable or disable profiling for a cmesh. If profiling is enabled, runtimes
 * and statistics are collected during cmesh_commit.
 * \param [in,out] cmesh          The cmesh to be updated.
 * \param [in]     set_profiling  If true, profiling will be enabled, if false
 *                                disabled.
 *
 * Profiling is disabled by default.
 * The cmesh must not be committed before calling this function.
 * \see t8_cmesh_print_profile
 */
void
t8_cmesh_set_profiling (t8_cmesh_t cmesh, int set_profiling);

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
int
t8_cmesh_is_equal (t8_cmesh_t cmesh_a, t8_cmesh_t cmesh_b);

/** Check whether a cmesh is empty on all processes.
 * \param [in]  cmesh           A committed cmesh.
 * \return                      True (non-zero) if and only if the cmesh has trees at all.
 */
int
t8_cmesh_is_empty (t8_cmesh_t cmesh);

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
 * \note It is illegal to broadcast a cmesh with a registered geometry (\ref t8_cmesh_register_geometry).
 *       All geometries must be registered after the broadcast (You can set tree attributes before bcast, though).
 */
t8_cmesh_t
t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm);

#if T8_ENABLE_METIS
/* TODO: document this. */
/* TODO: think about making this a pre-commit set_reorder function. */
void
t8_cmesh_reorder (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/* TODO: think about a sensible interface for a parmetis reordering. */
#endif

/** Register a geometry in the cmesh. The cmesh takes ownership of the geometry.
 * \param [in,out] cmesh        The cmesh.
 * \param [in]     geometry     The geometry to register.
 * 
 * If no geometry is registered and cmesh is modified from another cmesh then
 * the other cmesh's geometries are used.
 * \note If you need to use \ref t8_cmesh_bcast, then all geometries must be
 *       registered \a after the bcast operation, not before.
 */
void
t8_cmesh_register_geometry (t8_cmesh_t cmesh, t8_geometry_c *geometry);

/** Set the geometry for a tree, thus specify which geometry to use for this tree.
 * \param [in] cmesh     A non-committed cmesh.
 * \param [in] gtreeid   A global tree id in \a cmesh.
 * \param [in] geom      The geometry to use for this tree.
 * See also \ref t8_cmesh_get_tree_geometry
 */
void
t8_cmesh_set_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const t8_geometry_c *geom);

/** After allocating and adding properties to a cmesh, finish its construction.
 * TODO: this function is MPI collective.
 * \param [in,out] cmesh        Must be created with \ref t8_cmesh_init
 *                              (TODO: or bcast) and
 *                              specialized with t8_cmesh_set_* calls first (?).
 */
void
t8_cmesh_commit (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/* TODO: Document */
/* Currently, it is only legal to save cmeshes that use the linear geometry. */
int
t8_cmesh_save (t8_cmesh_t cmesh, const char *fileprefix);

/* TODO: Document */
t8_cmesh_t
t8_cmesh_load (const char *filename, sc_MPI_Comm comm);

/* TODO: Document */
/* procs_per_node is only relevant in mode==JUQUEEN.
 *  num_files = 1 => replicated cmesh is constructed */
t8_cmesh_t
t8_cmesh_load_and_distribute (const char *fileprefix, int num_files, sc_MPI_Comm comm, t8_load_mode_t mode,
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
int
t8_cmesh_comm_is_valid (t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Query whether a committed cmesh is partitioned or replicated.
 * \param [in] cmesh       A committed cmesh.
 * \return                 True if \a cmesh is partitioned.
 *                         False otherwise.
 * \a cmesh must be committed before calling this function.
 */
int
t8_cmesh_is_partitioned (t8_cmesh_t cmesh);

/** Get the dimension of a cmesh.
 * \param [in]  cmesh   The cmesh.
 * \a cmesh must be committed before calling this function.
 */
int
t8_cmesh_get_dimension (const t8_cmesh_t cmesh);

/** Return the global number of trees in a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of trees associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_gloidx_t
t8_cmesh_get_num_trees (t8_cmesh_t cmesh);

/** Return the number of local trees of a cmesh.
 *  If the cmesh is not partitioned this is equivalent to \ref t8_cmesh_get_num_trees.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of local trees of the cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t
t8_cmesh_get_num_local_trees (t8_cmesh_t cmesh);

/** Return the number of ghost trees of a cmesh.
 *  If the cmesh is not partitioned this is equivalent to \ref t8_cmesh_get_num_trees.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The number of ghost trees of the cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t
t8_cmesh_get_num_ghosts (t8_cmesh_t cmesh);

/** Return the global index of the first local tree of a cmesh.
 * If the cmesh is not partitioned this is always 0.
 * \param [in] cmesh       The cmesh to be considered.
 * \return                 The global id of the first local tree in cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_gloidx_t
t8_cmesh_get_first_treeid (t8_cmesh_t cmesh);

/** Get the geometry of a tree.
 * \param [in] cmesh   The cmesh.
 * \param [in] gtreeid The global tree id of the tree for which the geometry should be returned.
 * \return             The geometry of the tree.
 */
const t8_geometry_c *
t8_cmesh_get_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

/** Query whether a given t8_locidx_t belongs to a local tree of a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     An (possible) tree index.
 * \return                 True if \a ltreeid matches the range of local trees of \a cmesh.
 *                         False if not.
 * \a cmesh must be committed before calling this function.
 */
int
t8_cmesh_treeid_is_local_tree (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid);

/** Query whether a given t8_locidx_t belongs to a ghost of a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     An (possible) ghost index.
 * \return                 True if \a ltreeid matches the range of ghost trees of \a cmesh.
 *                         False if not.
 * \a cmesh must be committed before calling this function.
 */
int
t8_cmesh_treeid_is_ghost (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid);

/** Given a local tree id that belongs to a ghost, return the index of the ghost.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] ltreeid     The local id of a ghost, satisfying \ref t8_cmesh_treeid_is_ghost,
 *                         thus num_trees <= \a ltreeid < num_trees + num_ghosts
 * \return                 The index of the ghost within all ghosts, thus an index
 *                         0 <= index < num_ghosts
 * \a cmesh must be committed before calling this function.
 */
t8_locidx_t
t8_cmesh_ltreeid_to_ghostid (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid);

/* TODO: Replace this iterator with a new one that does not need the
 *        treeid to be part of the ctree struct */
/* TODO: should this and the next function be part of the interface? */
/** Return a pointer to the first local tree in a cmesh.
 * \param [in]     cmesh        The cmesh to be queried.
 * \return                      A pointer to the first local tree in \a cmesh.
 *                              If \a cmesh has no local trees, NULL is returned.
 * \a cmesh must be committed before calling this function.
 */
t8_ctree_t
t8_cmesh_get_first_tree (t8_cmesh_t cmesh);

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
t8_ctree_t
t8_cmesh_get_next_tree (t8_cmesh_t cmesh, t8_ctree_t tree);

/** Return a pointer to a given local tree.
 * \param [in]     cmesh        The cmesh to be queried.
 * \param [in]     ltree_id     The local id of the tree that is asked for.
 * \return                      A pointer to tree in \a cmesh with local
 *                              id \a ltree_id.
 * The cmesh must have at least \a ltree_id + 1 local trees when
 * calling this function.
 * \a cmesh must be committed before calling this function.
 */
t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_locidx_t ltree_id);

/** Return the eclass of a given local tree.
 * TODO: Should we refer to indices or consequently use ctree_t?
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    tree_id       The local id of the tree whose eclass will be returned.
 * \return                      The eclass of the given tree.
 * TODO: Call tree ids ltree_id or gtree_id etc. instead of tree_id.
 * \a cmesh must be committed before calling this function.
 */
t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_locidx_t ltree_id);

/** Query whether a face of a local tree or ghost is at the domain boundary.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    ltree_id      The local id of a tree.
 * \param [in]    face          The number of a face of the tree.
 * \return                      True if the face is at the domain boundary.
 *                              False otherwise.
 * \a cmesh must be committed before calling this function.
 */
int
t8_cmesh_tree_face_is_boundary (t8_cmesh_t cmesh, t8_locidx_t ltree_id, int face);

/** Return the eclass of a given local ghost.
 * TODO: Should we refer to indices or consequently use cghost_t?
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    ghost_id      The local id of the ghost whose eclass will be returned.
 *                              0 <= \a tree_id < cmesh.num_ghosts.
 * \return                      The eclass of the given ghost.
 * \a cmesh must be committed before calling this function.
 */
t8_eclass_t
t8_cmesh_get_ghost_class (t8_cmesh_t cmesh, t8_locidx_t lghost_id);

/** Return the global id of a given local tree or ghost.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    local_id      The local id of a tree or a ghost.
 *                              If \a local_id < cmesh.num_local_trees then it is
 *                              a tree, otherwise a ghost.
 * \return                      The global id of the tree/ghost.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_gloidx_t
t8_cmesh_get_global_id (t8_cmesh_t cmesh, t8_locidx_t local_id);

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
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_locidx_t
t8_cmesh_get_local_id (t8_cmesh_t cmesh, t8_gloidx_t global_id);

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
t8_locidx_t
t8_cmesh_get_face_neighbor (const t8_cmesh_t cmesh, const t8_locidx_t ltreeid, const int face, int *dual_face,
                            int *orientation);

/** Print the collected statistics from a cmesh profile.
 * \param [in]    cmesh         The cmesh.
 *
 * \a cmesh must be committed before calling this function.
 * \see t8_cmesh_set_profiling
 */
void
t8_cmesh_print_profile (t8_cmesh_t cmesh);

/** Return a pointer to the vertex coordinates of a tree.
 * \param [in]    cmesh         The cmesh.
 * \param [in]    ltreeid       The id of a local tree.
 * \return    If stored, a pointer to the vertex coordinates of \a tree.
 *            If no coordinates for this tree are found, NULL.
 */
double *
t8_cmesh_get_tree_vertices (t8_cmesh_t cmesh, t8_locidx_t ltreeid);

/** Return the attribute pointer of a tree.
 * \param [in]     cmesh        The cmesh.
 * \param [in]     package_id   The identifier of a valid software package. \see sc_package_register
 * \param [in]     key          A key used to identify the attribute under all
 *                              attributes of this tree with the same \a package_id.
 * \param [in]     tree_id      The local number of the tree.
 * \return         The attribute pointer of the tree \a ltree_id or NULL if the attribute is not found.
 * \note \a cmesh must be committed before calling this function.
 * \see t8_cmesh_set_attribute
 */
void *
t8_cmesh_get_attribute (const t8_cmesh_t cmesh, const int package_id, const int key, const t8_locidx_t ltree_id);

/** Return the attribute pointer of a tree for a gloidx_t array.
 * \param [in]     cmesh        The cmesh.
 * \param [in]     package_id   The identifier of a valid software package. \see sc_package_register
 * \param [in]     key          A key used to identify the attribute under all
 *                              attributes of this tree with the same \a package_id.
 * \param [in]     ltree_id     The local number of the tree.
 * \param [in]     data_count   The number of entries in the array that are requested. 
 *                              This must be smaller or equal to the \a data_count parameter
 *                              of the corresponding call to \ref t8_cmesh_set_attribute_gloidx_array
 * \return         The attribute pointer of the tree \a ltree_id or NULL if the attribute is not found.
 * \note \a cmesh must be committed before calling this function.
 * \note No check is performed whether the attribute actually stored \a data_count many entries since
 *       we do not store the number of data entries of the attribute array.
 *       You can keep track of the data count yourself by using another attribute.
 * \see t8_cmesh_set_attribute_gloidx_array
 */
t8_gloidx_t *
t8_cmesh_get_attribute_gloidx_array (const t8_cmesh_t cmesh, const int package_id, const int key,
                                     const t8_locidx_t ltree_id, const size_t data_count);

/** Return the shared memory array storing the partition table of
 * a partitioned cmesh.
 * \param [in]      cmesh       The cmesh.
 * \return                      The partition array.
 *                              NULL if the cmesh is not partitioned or
 *                              the partition array is not stored in \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
t8_shmem_array_t
t8_cmesh_get_partition_table (t8_cmesh_t cmesh);

/* TODO: remove get_ when there is no risk of confusion? Convention?
 *       Update: use get throughout for access functions that do not change the object.
 * */

/** Calculate the section of a uniform forest for the current rank.
 * \param [in]    cmesh         The cmesh to be considered.
 * \param [in]    level         The uniform refinement level to be created.
 * \param [in]    scheme            The element scheme for which to compute the bounds.
 * \param [out]   first_local_tree  The first tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_begin The global index of the first element belonging to the calling processor. Not computed if NULL.
 * \param [out]   last_local_tree  The last tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_end The global index of the first element that does not belonging to
 *                                  the calling processor anymore. Not computed if NULL.
 * \param [out]   first_tree_shared If not NULL, 1 or 0 is stored here depending on whether \a first_local_tree is the
 *                                 same as \a last_local_tree on the next process.
 * \a cmesh must be committed before calling this function. 
 */
void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, const int level, const t8_scheme_c *scheme, t8_gloidx_t *first_local_tree,
                         t8_gloidx_t *child_in_tree_begin, t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                         int8_t *first_tree_shared);

/**
 * Calculate the section of a uniform hybrid forest for the current rank. Need for hybrid meshes, especially 
 * meshes where not all elements refine into 1:2^dim manner. The section is calculated without assuming such refinement
 * and each process computes its number of elements on the given \var level communicates the number to other processes
 * and the correct section is computed based on this information. 
 * 
 * \param [in] cmesh        The cmesh to be considered.
 * \param [in] level        The uniform refinement level to be created.
 * \param [in] scheme       The element scheme for which to compute the bounds.
 * \param [out]   first_local_tree  The global index of the first tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_begin The global index of the first element belonging to the calling processor. Not computed if NULL.
 * \param [out]   last_local_tree  The global index of the last tree that contains elements belonging to the calling processor.
 * \param [out]   child_in_tree_end The global index of the first element that does not belonging to
 *                                  the calling processor anymore. Not computed if NULL.
 * \param [out]   first_tree_shared If not NULL, 1 or 0 is stored here depending on whether \a first_local_tree is the
 *                                 same as \a last_local_tree on the next process.
 * \param [in] comm         The communicator 
 */
void
t8_cmesh_uniform_bounds_for_irregular_refinement (const t8_cmesh_t cmesh, const int level, const t8_scheme_c *scheme,
                                                  t8_gloidx_t *first_local_tree, t8_gloidx_t *child_in_tree_begin,
                                                  t8_gloidx_t *last_local_tree, t8_gloidx_t *child_in_tree_end,
                                                  int8_t *first_tree_shared, sc_MPI_Comm comm);

/** Increase the reference counter of a cmesh.
 * \param [in,out] cmesh        On input, this cmesh must exist with positive
 *                              reference count.  It may be in any state.
 */
void
t8_cmesh_ref (t8_cmesh_t cmesh);

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
void
t8_cmesh_unref (t8_cmesh_t *pcmesh);

/** Verify that a coarse mesh has only one reference left and destroy it.
 * This function is preferred over \ref t8_cmesh_unref when it is known
 * that the last reference is to be deleted.
 * \param [in,out]  pcmesh      This cmesh must have a reference count of one.
 *                              It can be in any state (committed or not).
 *                              Then it effectively calls \ref t8_cmesh_unref.
 * \param [in]      comm        A mpi communicator that is valid with \a cmesh.
 */
void
t8_cmesh_destroy (t8_cmesh_t *pcmesh);

/** Compute y = ax + b on an array of doubles, interpreting
 * each 3 as one vector x 
 * \param[in]   coords_in         The incoming coordinates of the vectors
 * \param[out]  coords_out        The computed coordinates of the vectors
 * \param[in]   num_vertices      The number of vertices/vectors
 * \param[in]   alpha             Scaling factor for the vectors
 * \param[in]   b                 Translation of the vectors.*/

void
t8_cmesh_coords_axb (const double *coords_in, double *coords_out, int num_vertices, double alpha, const double b[3]);

/** Compute y = x + translate on an array of doubles, interpreting 
 * each 3 as one vector x
 * \param[in]   coords_in         The incoming coordinates of the vectors
 * \param[out]  coords_out        The computed coordinates of the vectors
 * \param[in]   num_vertices      The number of vertices/vectors
 * \param[in]   translate         Translation of the vectors.
 */
void
t8_cmesh_translate_coordinates (const double *coords_in, double *coords_out, const int num_vertices,
                                const double translate[3]);

/**TODO: Add proper documentation*/
void
t8_cmesh_new_translate_vertices_to_attributes (const t8_locidx_t *tvertices, const double *vertices,
                                               double *attr_vertices, const int num_vertices);

/**
 * \warning This function is only available in debug-modus and should only 
 * be used in debug-modus.
 * 
 * Prints the vertices of each tree of each process
 * 
 * \param[in] cmesh   Source-cmesh, which trees get printed.
 */
void
t8_cmesh_debug_print_trees (const t8_cmesh_t cmesh, sc_MPI_Comm comm);
T8_EXTERN_C_END ();

#endif /* !T8_CMESH_H */
