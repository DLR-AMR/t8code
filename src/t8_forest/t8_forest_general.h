/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_forest_general.h
 * We define the forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_GENERAL_H
#define T8_FOREST_GENERAL_H

#include <t8_cmesh/t8_cmesh.h>
#include <t8_element.h>
#include <t8_data/t8_containers.h>

/** Opaque pointer to a forest implementation. */
typedef struct t8_forest *t8_forest_t;

/** Opaque pointer to a tree implementation. */
typedef struct t8_tree *t8_tree_t;

/** This type controls, which neighbors count as ghost elements.
 * Currently, we support face-neighbors. Vertex and edge neighbors will eventually be added. */
typedef enum {
  T8_GHOST_NONE = 0, /**< Do not create ghost layer. */
  T8_GHOST_FACES,    /**< Consider all face (codimension 1) neighbors. */
  T8_GHOST_EDGES,    /**< Consider all edge (codimension 2) and face neighbors. */
  T8_GHOST_VERTICES  /**< Consider all vertex (codimension 3) and edge and face neighbors. */
} t8_ghost_type_t;

/** This typedef is needed as a helper construct to 
 * properly be able to define a function that returns
 * a pointer to a void fun(void) function. \see t8_forest_get_user_function.
 */
typedef void (*t8_generic_function_pointer) (void);
T8_EXTERN_C_BEGIN ();

/** Callback function prototype to replace one set of elements with another.
 *
 * This is used by the replace routine which can be called after adapt,
 * when the elements of an existing, valid
 * forest are changed. The callback allows the user to make changes to the elements
 * of the new forest that are either refined, coarsened or the same as elements in the old forest.
 *
 * \param [in] forest_old      The forest that is adapted
 * \param [in, out] forest_new The forest that is newly constructed from \a forest_old
 * \param [in] which_tree      The local tree containing \a first_outgoing and \a first_incoming
 * \param [in] tree_class      The eclass of the local tree containing \a first_outgoing and \a first_incoming
 * \param [in] scheme          The scheme of the forest
 * \param [in] refine          -1 if family in \a forest_old got coarsened, 0 if element
 *                             has not been touched, 1 if element got refined and -2 if
 *                             element got removed. See return of t8_forest_adapt_t.
 * \param [in] num_outgoing    The number of outgoing elements.
 * \param [in] first_outgoing  The tree local index of the first outgoing element.
 *                             0 <= first_outgoing < which_tree->num_elements
 * \param [in] num_incoming    The number of incoming elements.
 * \param [in] first_incoming  The tree local index of the first incoming element.
 *                             0 <= first_incom < new_which_tree->num_elements
 *
 * If an element is being refined, \a refine and \a num_outgoing will be 1 and 
 * \a num_incoming will be the number of children.
 * If a family is being coarsened, \a refine will be -1, \a num_outgoing will be 
 * the number of family members and \a num_incoming will be 1. 
 * If an element is being removed, \a refine and \a num_outgoing will be 1 and 
 * \a num_incoming will be 0. 
 * Else \a refine will be 0 and \a num_outgoing and \a num_incoming will both be 1.
 * \see t8_forest_iterate_replace
 */
typedef void (*t8_forest_replace_t) (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                     const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
                                     const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
                                     const t8_locidx_t first_incoming);

/** Callback function prototype to decide for refining and coarsening.
 * If \a is_family equals 1, the first \a num_elements in \a elements
 * form a family and we decide whether this family should be coarsened
 * or only the first element should be refined.
 * Otherwise \a is_family must equal zero and we consider the first entry
 * of the element array for refinement. 
 * Entries of the element array beyond the first \a num_elements are undefined.
 * \param [in] forest       The forest to which the new elements belong.
 * \param [in] forest_from  The forest that is adapted.
 * \param [in] which_tree   The local tree containing \a elements.
 * \param [in] tree_class   The eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in \a forest_from in the tree of the current element.
 * \param [in] scheme       The scheme of the forest.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *        -2 if the first entry in \a elements should be removed,
 *         0 else.
 */
/* TODO: Do we really need the forest argument? Since the forest is not committed yet it
 *       seems dangerous to expose to the user. */
typedef int (*t8_forest_adapt_t) (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                                  const t8_eclass_t tree_class, t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                                  const int is_family, const int num_elements, t8_element_t *elements[]);

/** Create a new forest with reference count one.
 * This forest needs to be specialized with the t8_forest_set_* calls.
 * Currently it is mandatory to either call the functions \see t8_forest_set_mpicomm, 
 * \ref t8_forest_set_cmesh, and \ref t8_forest_set_scheme,
 * or to call one of \ref t8_forest_set_copy, \ref t8_forest_set_adapt, or
 * \ref t8_forest_set_partition.  It is illegal to mix these calls, or to
 * call more than one of the three latter functions
 * Then it needs to be set up with \ref t8_forest_commit.
 * \param [in,out] pforest      On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new forest.
 */
void
t8_forest_init (t8_forest_t *pforest);

/** Check whether a forest is not NULL, initialized and not committed.
 * In addition, it asserts that the forest is consistent as much as possible.
 * \param [in] forest           This forest is examined.  May be NULL.
 * \return                      True if forest is not NULL,
 *                              \ref t8_forest_init has been called on it,
 *                              but not \ref t8_forest_commit.
 *                              False otherwise.
 */
int
t8_forest_is_initialized (t8_forest_t forest);

/** Check whether a forest is not NULL, initialized and committed.
 * In addition, it asserts that the forest is consistent as much as possible.
 * \param [in] forest           This forest is examined.  May be NULL.
 * \return                      True if forest is not NULL and
 *                              \ref t8_forest_init has been called on it
 *                              as well as \ref t8_forest_commit.
 *                              False otherwise.
 */
int
t8_forest_is_committed (t8_forest_t forest);

/** Check whether the forest has local overlapping elements.
 * \param [in] forest   The forest to consider.
 * \return              True if \a forest has no elements which are inside each other.
 * \note This function is collective, but only checks local overlapping on each process.
 * \see t8_forest_partition_test_boundary_element if you also want to test for 
 * global overlap across the process boundaries.
 */
int
t8_forest_no_overlap (t8_forest_t forest);

/** Check whether two committed forests have the same local elements.
 * \param [in] forest_a The first forest.
 * \param [in] forest_b The second forest.
 * \return              True if \a forest_a and \a forest_b do have the same
 *                      number of local trees and each local tree has the same
 *                      elements, that is \ref t8_element_is_equal returns true
 *                      for each pair of elements of \a forest_a and \a forest_b.
 * \note This function is not collective. It only returns the state on the current
 * rank.
 */
int
t8_forest_is_equal (t8_forest_t forest_a, t8_forest_t forest_b);

/** Set the cmesh associated to a forest.
 * By default, the forest takes ownership of the cmesh such that it will be
 * destroyed when the forest is destroyed.  To keep ownership of the cmesh,
 * call \ref t8_cmesh_ref before passing it to \ref t8_forest_set_cmesh.
 * This means that it is ILLEGAL to continue using cmesh or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest       The forest whose cmesh variable will be set.
 * \param [in]     cmesh        The cmesh to be set.  We take ownership.
 *                              This can be prevented by referencing \b cmesh.
 * \param [in]     comm         The MPI communicator.
 */
void
t8_forest_set_cmesh (t8_forest_t forest, t8_cmesh_t cmesh, sc_MPI_Comm comm);

/** Set the element scheme associated to a forest.
 * By default, the forest takes ownership of the scheme such that it will be
 * destroyed when the forest is destroyed.  To keep ownership of the scheme, call
 * \ref t8_scheme_ref before passing it to \ref t8_forest_set_scheme.
 * This means that it is ILLEGAL to continue using scheme or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest       The forest whose scheme variable will be set.
 * \param [in]     scheme       The scheme to be set.  We take ownership.
 *                              This can be prevented by referencing \b scheme.
 */
void
t8_forest_set_scheme (t8_forest_t forest, const t8_scheme_c *scheme);

/** Set the initial refinement level to be used when \b forest is committed.
 * \param [in,out] forest      The forest whose level will be set.
 * \param [in]     level       The initial refinement level of \b forest, when
 *                             it is committed.
 * \note This setting cannot be combined with any of the derived forest methods
 * (\ref t8_forest_set_copy, \ref t8_forest_set_adapt, \ref t8_forest_set_partition,
 * and \ref t8_forest_set_balance) and overwrites any of these settings.
 * If this function is used, then the forest is created from scratch as a uniform
 * refinement of the specified cmesh (\ref t8_forest_set_cmesh, \ref t8_forest_set_scheme).
 */
void
t8_forest_set_level (t8_forest_t forest, int level);

/** Set a forest as source for copying on committing.
 * By default, the forest takes ownership of the source \b from such that it will
 * be destroyed on calling \ref t8_forest_commit.  To keep ownership of \b
 * from, call \ref t8_forest_ref before passing it into this function.
 * This means that it is ILLEGAL to continue using \b from or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest     The forest.
 * \param [in]     from       A second forest from which \a forest will be copied
 *                            in \ref t8_forest_commit.
 * \note This setting cannot be combined with \ref t8_forest_set_adapt,
 * \ref t8_forest_set_partition, or \ref t8_forest_set_balance and overwrites these
 * settings.
 */
void
t8_forest_set_copy (t8_forest_t forest, const t8_forest_t from);

/** Set a source forest with an adapt function to be adapted on committing.
 * By default, the forest takes ownership of the source \b set_from such that it
 * will be destroyed on calling \ref t8_forest_commit. To keep ownership of \b
 * set_from, call \ref t8_forest_ref before passing it into this function.
 * This means that it is ILLEGAL to continue using \b set_from or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest   The forest
 * \param [in] set_from     The source forest from which \b forest will be adapted.
 *                          We take ownership. This can be prevented by
 *                          referencing \b set_from.
 *                          If NULL, a previously (or later) set forest will
 *                          be taken (\ref t8_forest_set_partition, \ref t8_forest_set_balance).
 * \param [in] adapt_fn     The adapt function used on committing.
 * \param [in] recursive    A flag specifying whether adaptation is to be done recursively
 *                          or not. If the value is zero, adaptation is not recursive
 *                          and it is recursive otherwise.
 * \note This setting can be combined with \ref t8_forest_set_partition and \ref
 * t8_forest_set_balance. The order in which these operations are executed is always
 * 1) Adapt 2) Partition 3) Balance.
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
/* TODO: make recursive flag to int specifying the number of recursions? */
void
t8_forest_set_adapt (t8_forest_t forest, const t8_forest_t set_from, t8_forest_adapt_t adapt_fn, int recursive);

/** Set the user data of a forest. This can i.e. be used to pass user defined
 * arguments to the adapt routine.
 * \param [in,out] forest   The forest
 * \param [in]     data     A pointer to user data. t8code will never touch the data.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_get_user_data
 */
void
t8_forest_set_user_data (t8_forest_t forest, void *data);

/** Return the user data pointer associated with a forest.
 * \param [in]     forest   The forest.
 * \return                  The user data pointer of \a forest.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_set_user_data
 */
void *
t8_forest_get_user_data (const t8_forest_t forest);

/** Set the user function pointer of a forest. This can i.e. be used to pass user defined
 * functions to the adapt routine.
 * \param [in,out] forest   The forest
 * \param [in]     function A pointer to a user defined function. t8code will never touch the function.
 * The forest does not need be committed before calling this function.
 * \note \a function can be an arbitrary function with return value and parameters of
 * your choice. When accessing it with \ref t8_forest_get_user_function you should cast
 * it into the proper type.
 * \see t8_forest_get_user_function
 */
void
t8_forest_set_user_function (t8_forest_t forest, t8_generic_function_pointer function);

/** Return the user function pointer associated with a forest.
 * \param [in]     forest   The forest.
 * \return                  The user function pointer of \a forest.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_set_user_function
 */
t8_generic_function_pointer
t8_forest_get_user_function (const t8_forest_t forest);

/** Set a source forest to be partitioned during commit.
 * The partitioning is done according to the SFC and each rank is assigned
 * the same (maybe +1) number of elements.
 * \param [in, out] forest  The forest.
 * \param [in]      set_from A second forest that should be partitioned.
 *                          We take ownership. This can be prevented by
 *                          referencing \b set_from.
 *                          If NULL, a previously (or later) set forest will
 *                          be taken (\ref t8_forest_set_adapt, \ref t8_forest_set_balance).
 * \param [in]      set_for_coarsening CURRENTLY DISABLED. If true, then the partitions
 *                          are choose such that coarsening an element once is a process local
 *                          operation.
 * \note This setting can be combined with \ref t8_forest_set_adapt and \ref
 * t8_forest_set_balance. The order in which these operations are executed is always
 * 1) Adapt 2) Partition 3) Balance.
 * If \ref t8_forest_set_balance is called with the \a no_repartition parameter set as
 * false, it is not necessary to call \ref t8_forest_set_partition additionally.
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
void
t8_forest_set_partition (t8_forest_t forest, const t8_forest_t set_from, int set_for_coarsening);

/** Set a source forest to be balanced during commit.
 * A forest is said to be balanced if each element has face neighbors of level
 * at most +1 or -1 of the element's level.
 * \param [in, out] forest  The forest.
 * \param [in]      set_from A second forest that should be balanced.
 *                          We take ownership. This can be prevented by
 *                          referencing \b set_from.
 *                          If NULL, a previously (or later) set forest will
 *                          be taken (\ref t8_forest_set_adapt, \ref t8_forest_set_partition)
 * \param [in]      no_repartition Balance constructs several intermediate forest that
 *                          are refined from each other. In order to maintain a balanced load
 *                          these forest are repartitioned in each round and the resulting
 *                          forest is load-balanced per default.
 *                          If this behaviour is not desired, \a no_repartition should be
 *                          set to true.
 *                          If \a no_repartition is false, an additional call of \ref t8_forest_set_partition is not
 *                          necessary.
 * \note This setting can be combined with \ref t8_forest_set_adapt and \ref
 * t8_forest_set_partition. The order in which these operations are executed is always
 * 1) Adapt 2) Partition 3) Balance.
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
void
t8_forest_set_balance (t8_forest_t forest, const t8_forest_t set_from, int no_repartition);

/** Enable or disable the creation of a layer of ghost elements.
 * On default no ghosts are created.
 * \param [in]      forest    The forest.
 * \param [in]      do_ghost  If non-zero a ghost layer will be created.
 * \param [in]      ghost_type Controls which neighbors count as ghost elements,
 *                             currently only T8_GHOST_FACES is supported. This value
 *                             is ignored if \a do_ghost = 0.
 */
void
t8_forest_set_ghost (t8_forest_t forest, int do_ghost, t8_ghost_type_t ghost_type);

/** Like \ref t8_forest_set_ghost but with the additional options to change the
 * ghost algorithm. This is used for debugging and timing the algorithm.
 * An application should almost always use \ref t8_forest_set_ghost.
 * \param [in]      forest        The forest.
 * \param [in]      do_ghost      If non-zero a ghost layer will be created.
 * \param [in]      ghost_type    Controls which neighbors count as ghost elements,
 *                                currently only T8_GHOST_FACES is supported. This value
 *                                is ignored if \a do_ghost = 0.
 * \param [in]      ghost_version If 1, the iterative ghost algorithm for balanced forests is used.
 *                                If 2, the iterative algorithm for unbalanced forests.
 *                                If 3, the top-down search algorithm for unbalanced forests.
 * \see t8_forest_set_ghost
 */
void
t8_forest_set_ghost_ext (t8_forest_t forest, int do_ghost, t8_ghost_type_t ghost_type, int ghost_version);

/**
 *  Use assertions and document that the forest_set (..., from) and
 *  set_load are mutually exclusive. 
 * 
 *  TODO: Unused function -> remove?
 */
void
t8_forest_set_load (t8_forest_t forest, const char *filename);

/** Compute the global number of leaf elements in a forest as the sum
 *  of the local leaf element counts.
 *  \param [in] forest    The forest.
 */
void
t8_forest_comm_global_num_leaf_elements (t8_forest_t forest);

/** After allocating and adding properties to a forest, commit the changes.
 * This call sets up the internal state of the forest.
 * \param [in,out] forest       Must be created with \ref t8_forest_init and
 *                              specialized with t8_forest_set_* calls first.
 */
void
t8_forest_commit (t8_forest_t forest);

/** Return the maximum allowed refinement level for any element in a forest.
 * \param [in]  forest    A forest.
 * \return                The maximum level of refinement that is allowed for
 *                        an element in this forest. It is guaranteed that any tree
 *                        in \a forest can be refined this many times and it is not
 *                        allowed to refine further.
 * \a forest must be committed before calling this function.
 * For forest with a single element class (non-hybrid) maxlevel is the maximum
 * refinement level of this element class, whilst for hybrid forests the maxlevel is
 * the minimum of all maxlevels of the element classes in this forest.
 */
int
t8_forest_get_maxlevel (const t8_forest_t forest);

/** Return the number of process local leaf elements in the forest.
 * \param [in]  forest    A forest.
 * \return                The number of leaf elements on this process in \a forest.
 * \a forest must be committed before calling this function.
 */
t8_locidx_t
t8_forest_get_local_num_leaf_elements (const t8_forest_t forest);

/** Return the number of global leaf elements in the forest.
 * \param [in]  forest    A forest.
 * \return                The number of leaf elements (summed over all processes) in \a forest.
 * \a forest must be committed before calling this function.
 */
t8_gloidx_t
t8_forest_get_global_num_leaf_elements (const t8_forest_t forest);

/** Return the number of ghost elements of a forest.
 * \param [in]      forest      The forest.
 * \return                      The number of ghost elements stored in the ghost
 *                              structure of \a forest. 0 if no ghosts were constructed.
 *                              \see t8_forest_set_ghost
 * \a forest must be committed before calling this function.
 */
t8_locidx_t
t8_forest_get_num_ghosts (const t8_forest_t forest);

/** Return the element class of a forest local tree.
 * \param [in] forest    The forest.
 * \param [in] ltreeid   The local id of a tree in \a forest.
 * \return  The element class of the tree \a ltreeid.
 * \a forest must be committed before calling this function.
 */
t8_eclass_t
t8_forest_get_eclass (const t8_forest_t forest, const t8_locidx_t ltreeid);

/**
 * Check whether a given tree id belongs to a local tree in a forest.
 * 
 * \param [in]    forest The forest.
 * \param [in]    local_tree A tree id.
 * \return True if and only if the id \a local_tree belongs to a local tree of \a forest.
 * \a forest must be committed before calling this function.
 */
int
t8_forest_tree_is_local (const t8_forest_t forest, const t8_locidx_t local_tree);

/** Given a global tree id compute the forest local id of this tree.
 * If the tree is a local tree, then the local id is between 0 and the number
 * of local trees. If the tree is not a local tree, a negative number is returned.
 * \param [in]      forest The forest.
 * \param [in]      gtreeid The global id of a tree.
 * \return                 The tree's local id in \a forest, if it is a local tree.
 *                         A negative number if not. Ghosts trees are not considered 
 *                         as local. 
 * \see t8_forest_get_local_or_ghost_id for ghost trees.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_locidx_t
t8_forest_get_local_id (const t8_forest_t forest, const t8_gloidx_t gtreeid);

/** Given a global tree id compute the forest local id of this tree.
 * If the tree is a local tree, then the local id is between 0 and the number
 * of local trees. If the tree is a ghost, then the local id is between num_local_trees and
 * num_local_trees + num_ghost_trees.
 * If the tree is neither a local tree nor a ghost tree, a negative number is returned.
 * \param [in]      forest The forest.
 * \param [in]      gtreeid The global id of a tree.
 * \return                 The tree's local id in \a forest, if it is a local tree.
 *                         num_local_trees + the ghosts id, if it is a ghost tree.
 *                         A negative number if not.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing
 */
t8_locidx_t
t8_forest_get_local_or_ghost_id (const t8_forest_t forest, const t8_gloidx_t gtreeid);

/** Given the local id of a tree in a forest, compute the tree's local id in the associated cmesh.
 * \param [in] forest    The forest.
 * \param [in] ltreeid   The local id of a tree or ghost in the forest.
 * \return  The local id of the tree in the cmesh associated with the forest.
 * \a forest must be committed before calling this function.
 * \note For forest local trees, this is the inverse function of \ref t8_forest_cmesh_ltreeid_to_ltreeid.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_locidx_t
t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest, t8_locidx_t ltreeid);

/** Given the local id of a tree in the coarse mesh of a forest, compute the tree's local id in the forest.
 * \param [in] forest    The forest.
 * \param [in] lctreeid  The local id of a tree in the coarse mesh of \a forest.
 * \return  The local id of the tree in the forest. -1 if the tree is not forest local.
 * \a forest must be committed before calling this function.
 * \note For forest local trees, this is the inverse function of \ref t8_forest_ltreeid_to_cmesh_ltreeid.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_locidx_t
t8_forest_cmesh_ltreeid_to_ltreeid (t8_forest_t forest, t8_locidx_t lctreeid);

/** Given the local id of a tree in a forest, return the coarse tree of the cmesh that corresponds to this tree.
 * \param [in] forest     The forest.
 * \param [in] ltreeid    The local id of a tree in the forest.
 * \return                The coarse tree that matches the forest tree with local id \a ltreeid.
 */
t8_ctree_t
t8_forest_get_coarse_tree (t8_forest_t forest, t8_locidx_t ltreeid);

/**
 * Query whether a given element is a leaf in a forest.
 * 
 * \param [in]  forest    The forest.
 * \param [in]  element   An element of a local tree in \a forest.
 * \param [in]  local_tree A local tree id of \a forest.
 * \return True (non-zero) if and only if \a element is a leaf in \a local_tree of \a forest.
 * \note This does not query for ghost leaves.
 * \note \a forest must be committed before calling this function.
 */
int
t8_forest_element_is_leaf (const t8_forest_t forest, const t8_element_t *element, const t8_locidx_t local_tree);

/**
 * Query whether a given element or a ghost is a leaf of a local or ghost tree in a forest.
 * 
 * \param [in]  forest    The forest.
 * \param [in]  element   An element of a local tree in \a forest.
 * \param [in]  local_tree A local tree id of \a forest or a ghost tree id
 * \param [in]  check_ghost If true \a element is interpreted as a ghost element and
 *                         \a local_tree as the id of a ghost tree (0 <= \a local_tree < num_ghost_trees).
 *                         If false \a element is interpreted as an element and \a local_tree as
 *                         the id of a local tree (0 <= \a local_tree < num_local_trees).
 * \return True (non-zero) if and only if \a element is a leaf (or ghost) in \a local_tree of \a forest.
 * \note \a forest must be committed before calling this function.
 * \ref t8_forest_element_is_leaf
 * \ref t8_forest_element_is_ghost
 */
int
t8_forest_element_is_leaf_or_ghost (const t8_forest_t forest, const t8_element_t *element, const t8_locidx_t local_tree,
                                    const int check_ghost);

/** Compute the leaf face orientation at given face in a forest.
 * \param [in]    forest  The forest. Must have a valid ghost layer.
 * \param [in]    ltreeid A local tree id.
 * \param [in]    scheme      The eclass scheme of the element.
 * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
 * \param [in]    face    The index of the face across which the face neighbors
 *                        are searched.
 * \return                Face orientation encoded as integer.
 *
 * For more information about the encoding of face orientation refer to \ref t8_cmesh_get_face_neighbor.
 */
int
t8_forest_leaf_face_orientation (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_scheme_c *scheme,
                                 const t8_element_t *leaf, const int face);

/** Compute the leaf face neighbors of a forest leaf element or ghost leaf.
 * \param [in]    forest  The forest.
 * \param [in]    ltreeid A local tree id (could also be a ghost tree). 0 <= \a ltreeid < num_local trees+num_ghost_trees
 * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
 * \param [out]   pneighbor_leaves Unallocated on input. On output the neighbor
 *                        leaves are stored here.
 * \param [in]    face    The index of the face across which the face neighbors
 *                        are searched.
 * \param [out]   dual_faces On output the face id's of the neighboring elements' faces.
 * \param [out]   num_neighbors On output the number of neighbor leaves.
 * \param [out]   pelement_indices Unallocated on input. On output the element indices
 *                        of the neighbor leaves are stored here.
 *                        0, 1, ... num_local_el - 1 for local leaves and
 *                        num_local_el , ... , num_local_el + num_ghosts - 1 for ghosts.
 * \param [out]   pneigh_eclass On output the eclass of the neighbor elements.
 * \note If there are no face neighbors, then *pneighbor_leaves = NULL, num_neighbors = 0,
 * and *pelement_indices = NULL on output.
 * \note \a forest must be committed before calling this function.
 * \note If \a forest does not have a ghost layer then leaf elements at the process boundaries have 0 neighbors. (The function output for leaf elements then depends on the parallel partition.)
 * \note Important! This routine allocates memory which must be freed. Do it like this:
 *
 *   if (num_neighbors > 0) {
 *     T8_FREE (pneighbor_leaves);
 *     T8_FREE (pelement_indices);
 *     T8_FREE (dual_faces);
 *   }
 *
 */
void
t8_forest_leaf_face_neighbors (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *leaf,
                               const t8_element_t **pneighbor_leaves[], const int face, int *dual_faces[],
                               int *num_neighbors, t8_locidx_t **pelement_indices, t8_eclass_t *pneigh_eclass);

/** Like \ref t8_forest_leaf_face_neighbors but also provides information about the global neighbors and the orientation. 
 * \param [in]    forest  The forest. Must have a valid ghost layer.
 * \param [in]    ltreeid A local tree id (could also be a ghost tree). 0 <= \a ltreeid < num_local trees+num_ghost_trees
 * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
 * \param [out]   pneighbor_leaves Unallocated on input. On output the neighbor
 *                        leaves are stored here.
 * \param [in]    face    The index of the face across which the face neighbors
 *                        are searched.
 * \param [out]   dual_faces On output the face id's of the neighboring elements' faces.
 * \param [out]   num_neighbors On output the number of neighbor leaves.
 * \param [out]   pelement_indices Unallocated on input. On output the element indices
 *                        of the neighbor leaves are stored here.
 *                        0, 1, ... num_local_el - 1 for local leaves and
 *                        num_local_el , ... , num_local_el + num_ghosts - 1 for ghosts.
 * \param [out]   pneigh_eclass On output the eclass of the neighbor elements.
 * \param [out]   gneigh_tree  The global tree IDs of the neighbor trees.
 * \param [out]   orientation  If not NULL on input, the face orientation is computed and stored here. 
 *                                         Thus, if the face connection is an inter-tree connection the orientation of the tree-to-tree connection is stored. 
 *                                         Otherwise, the value 0 is stored.
 * All other parameters and behavior are identical to \ref t8_forest_leaf_face_neighbors.
 * \note If there are no face neighbors, then *pneighbor_leaves = NULL, num_neighbors = 0,
 * and *pelement_indices = NULL on output.
 * \note \a forest must be committed before calling this function.
 *
 * \note Important! This routine allocates memory which must be freed. Do it like this:
 *
 *   if (num_neighbors > 0) {
 *     T8_FREE (pneighbor_leaves);
 *     T8_FREE (pelement_indices);
 *     T8_FREE (dual_faces);
 *   }
 *
 */
void
t8_forest_leaf_face_neighbors_ext (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                   const t8_element_t *leaf_or_ghost, const t8_element_t **pneighbor_leaves[],
                                   const int face, int *dual_faces[], int *num_neighbors,
                                   t8_locidx_t **pelement_indices, t8_eclass_t *pneigh_eclass, t8_gloidx_t *gneigh_tree,
                                   int *orientation);

/** Given a leaf element or ghost index in "all local elements + ghosts" enumeration
 * compute the index of the face neighbor of the element - provided that only one or no
 * face neighbors exists.
 * HANDLE WITH CARE. DO NOT CALL IF THE FOREST IS ADAPTED.
 * 
 * \param[in] forest        The forest. Must be committed.
 * \param[in] element_index Index of an element in \a forest. Must have only one or no facen neighbors across the given face.
 *                          0 <= \a element_index < num_local_elements + num_ghosts
 * \param[in] face_index    Index of a face of \a element.
 * \param[in] global_treeid Global index of the tree that contains \a element.
 * \param[out] dual_face    Return value, the dual_face index of the face neighbor.
 * \return The index of the face neighbor leaf (local element or ghost).
 * \note Do not call if you are unsure about the number of face neighbors. In particular if the forest is adapted and not uniform.
 */
t8_locidx_t
t8_forest_same_level_leaf_face_neighbor_index (const t8_forest_t forest, const t8_locidx_t element_index,
                                               const int face_index, const t8_gloidx_t global_treeid, int *dual_face);

/** Exchange ghost information of user defined element data.
 * \param [in] forest       The forest. Must be committed.
 * \param [in] element_data An array of length num_local_elements + num_ghosts
 *                         storing one value for each local element and ghost in \a forest.
 *                         After calling this function the entries for the ghost elements
 *                         are update with the entries in the \a element_data array of
 *                         the corresponding owning process.
 * \note This function is collective and hence must be called by all processes in the forest's
 *       MPI Communicator.
 */
/* TODO: In \ref t8_forest_ghost_cxx we already implemented a begin and end function
 *       that allow for overlapping communication and computation. We will make them
 *       available in this interface in the future. */
void
t8_forest_ghost_exchange_data (t8_forest_t forest, sc_array_t *element_data);

/** Print the ghost structure of a forest. Only used for debugging. */
void
t8_forest_ghost_print (t8_forest_t forest);

/** Change the cmesh associated to a forest to a partitioned cmesh that
 * is partitioned according to the tree distribution in the forest.
 * \param [in,out]   forest The forest.
 * \param [in]       comm   The MPI communicator that is used to partition
 *                          and commit the cmesh.
 * \param [in]       set_profiling If true, profiling for the new cmesh
 *                          will be enabled. \see t8_cmesh_set_profiling, \see t8_cmesh_print_profile
 *  \see t8_cmesh.h
 */
void
t8_forest_partition_cmesh (t8_forest_t forest, sc_MPI_Comm comm, int set_profiling);

/** Return the mpi communicator associated to a forest.
 * \param [in]      forest      The forest.
 * \return                      The mpi communicator of \a forest.
 * \a forest must be committed before calling this function.
 */
sc_MPI_Comm
t8_forest_get_mpicomm (const t8_forest_t forest);

/** Return the global id of the first local tree of a forest.
 * \param [in]      forest      The forest.
 * \return                      The global id of the first local tree in \a forest.
 */
t8_gloidx_t
t8_forest_get_first_local_tree_id (const t8_forest_t forest);

/** Return the number of local trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of local trees of that forest.
 */
t8_locidx_t
t8_forest_get_num_local_trees (const t8_forest_t forest);

/** Return the number of ghost trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of ghost trees of that forest.
 */
t8_locidx_t
t8_forest_get_num_ghost_trees (const t8_forest_t forest);

/** Return the number of global trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of global trees of that forest.
 */
t8_gloidx_t
t8_forest_get_num_global_trees (const t8_forest_t forest);

/** Return the global id of a local tree or a ghost tree.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     An id 0 <= \a ltreeid < num_local_trees + num_ghosts
 *                              specifying a local tree or ghost tree.
 * \return          The global id corresponding to the tree with local id \a ltreeid.
 * \a forest must be committed before calling this function.
 * \see https://github.com/DLR-AMR/t8code/wiki/Tree-indexing for more details about tree indexing.
 */
t8_gloidx_t
t8_forest_global_tree_id (const t8_forest_t forest, const t8_locidx_t ltreeid);

/** Return a pointer to a tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The local id of the tree.
 * \return                      A pointer to the tree with local id \a ltree_id.
 * \a forest must be committed before calling this function.
 */
t8_tree_t
t8_forest_get_tree (const t8_forest_t forest, const t8_locidx_t ltree_id);

/** Return a pointer to the vertex coordinates of a tree.
 * \param [in]    forest        The forest.
 * \param [in]    ltreeid       The id of a local tree.
 * \return    If stored, a pointer to the vertex coordinates of \a tree.
 *            If no coordinates for this tree are found, NULL.
 */
double *
t8_forest_get_tree_vertices (t8_forest_t forest, t8_locidx_t ltreeid);

/** Return the array of leaf elements of a local tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The local id of a local tree of \a forest.
 * \return                      An array of t8_element_t * storing all leaf elements
 *                              of this tree.
 */
t8_element_array_t *
t8_forest_tree_get_leaf_elements (const t8_forest_t forest, const t8_locidx_t ltree_id);

/** Return a cmesh associated to a forest.
 * \param [in]      forest      The forest.
 * \return          The cmesh associated to the forest.
 */
t8_cmesh_t
t8_forest_get_cmesh (t8_forest_t forest);

/** Return a leaf element of the forest.
 * \param [in]      forest      The forest.
 * \param [in]      lelement_id The local id of a leaf element in \a forest.
 * \param [out]     ltreeid     If not NULL, on output the local tree id of the tree in which the
 *                              leaf element lies in.
 * \return          A pointer to the leaf element. NULL if this element does not exist. Ghost elements are
 *                  not considered as local.
 * \see t8_forest_ghost_get_leaf_element to access ghost leaf elements.
 * \note This function performs a binary search. For constant access, use \ref t8_forest_get_leaf_element_in_tree
 * \a forest must be committed before calling this function.
 */
t8_element_t *
t8_forest_get_leaf_element (t8_forest_t forest, t8_locidx_t lelement_id, t8_locidx_t *ltreeid);

/** Return a leaf element of a local tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     An id of a local tree in the forest. Ghost trees are not considered local.
 * \param [in]      leid_in_tree The index of a leaf element in the tree.
 * \return          A pointer to the leaf element.
 * \see t8_forest_ghost_get_leaf_element_in_tree to access ghost leaf elements.
 * \note If the tree id is know, this function should be preferred over \ref t8_forest_get_leaf_element.
 * \a forest must be committed before calling this function.
 */
const t8_element_t *
t8_forest_get_leaf_element_in_tree (t8_forest_t forest, t8_locidx_t ltreeid, t8_locidx_t leid_in_tree);

/** Return the number of leaf elements of a tree.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     A local id of a tree.
 * \return                      The number of leaf elements in the local tree \a ltreeid.
 */
t8_locidx_t
t8_forest_get_tree_num_leaf_elements (t8_forest_t forest, t8_locidx_t ltreeid);

/** Return the element offset of a local tree, that is the number of leaf elements
 * in all trees with smaller local treeid.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     A local id of a tree.
 * \return                      The number of leaf elements on all local tree with
 *                              id < \a ltreeid.
 * \note \a forest must be committed before calling this function.
 */
t8_locidx_t
t8_forest_get_tree_element_offset (const t8_forest_t forest, const t8_locidx_t ltreeid);

/** Return the number of leaf elements of a tree.
 * \param [in]      tree       A tree in a forest.
 * \return                     The number of leaf elements of that tree.
 */
t8_locidx_t
t8_forest_get_tree_leaf_element_count (t8_tree_t tree);

/** Return the eclass of a tree in a forest.
 * \param [in]      forest    The forest.
 * \param [in]      ltreeid   The local id of a tree (local or ghost) in \a forest.
 * \return                    The element class of the tree with local id \a ltreeid.
 */
t8_eclass_t
t8_forest_get_tree_class (const t8_forest_t forest, const t8_locidx_t ltreeid);

/** Compute the global index of the first local leaf element of a forest.
 * This function is collective.
 * \param [in]     forest       A committed forest, whose first leaf element's index is computed.
 * \return         The global index of \a forest's first local leaf element.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 */
t8_gloidx_t
t8_forest_get_first_local_leaf_element_id (t8_forest_t forest);

/** Return the element scheme associated to a forest.
 * \param [in]      forest      A committed forest.
 * \return          The element scheme of the forest.
 * \see t8_forest_set_scheme
 */
const t8_scheme_c *
t8_forest_get_scheme (const t8_forest_t forest);

/** Return the eclass of the tree in which a face neighbor of a given element or ghost
 * lies.
 * \param [in]      forest      A committed forest.
 * \param [in]      ltreeid     The local tree or ghost tree in which the element lies. 0 <= \a ltreeid < num_local_trees + num_ghost_trees
 * \param [in]      elem        An element or ghost in the tree \a ltreeid.
 * \param [in]      face        A face number of \a elem.
 * \return                      The eclass of the local tree or ghost tree that
 *                              is face neighbor of \a elem across \a face.
 *                              T8_ECLASS_INVALID if no neighbor exists.
 */
t8_eclass_t
t8_forest_element_neighbor_eclass (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *elem,
                                   const int face);

/** Construct the face neighbor of an element, possibly across tree boundaries.
 * Returns the global tree-id of the tree in which the neighbor element lies in.
 *
 * \param [in] forest       The forest.
 * \param [in] ltreeid      The local tree in which the element lies.
 * \param [in] elem         The element to be considered.
 * \param [in,out] neigh    On input an allocated element of the scheme of the
 *                          face_neighbors eclass.
 *                          On output, this element's data is filled with the
 *                          data of the face neighbor. If the neighbor does not exist
 *                          the data could be modified arbitrarily.
 * \param [in] neigh_eclass The eclass of \a neigh.
 * \param [in] face         The number of the face along which the neighbor should be constructed.
 * \param [out] neigh_face  The number of the face viewed from perspective of \a neigh.
 * \return The global tree-id of the tree in which \a neigh is in.
 *        -1 if there exists no neighbor across that face. Domain boundary.
 *        -2 if the neighbor is not in a local tree or ghost tree. Process/Ghost boundary.
 */
t8_gloidx_t
t8_forest_element_face_neighbor (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem, t8_element_t *neigh,
                                 const t8_eclass_t neigh_eclass, int face, int *neigh_face);

/**
 * TODO: Can be removed since it is unused.
 * 
 * \param[in] forest The forest.
 */
void
t8_forest_iterate (t8_forest_t forest);

/** Query whether a batch of points lies inside an element. For bilinearly interpolated elements.
 * \note For 2D quadrilateral elements this function is only an approximation. It is correct
 *  if the four vertices lie in the same plane, but it may produce only approximate results if 
 *  the vertices do not lie in the same plane.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     The forest local id of the tree in which the element is.
 * \param [in]      element     The element.
 * \param [in]      points      3-dimensional coordinates of the points to check
 * \param [in]      num_points  The number of points to check
 * \param [in, out] is_inside   An array of length \a num_points, filled with 0/1 on output. True (non-zero) if a \a point 
 *                              lies within an \a element, false otherwise. The return value is also true if the point 
 *                              lies on the element boundary. Thus, this function may return true for different leaf 
 *                              elements, if they are neighbors and the point lies on the common boundary.
 * \param [in]      tolerance   Tolerance that we allow the point to not exactly match the element.
 *                              If this value is larger we detect more points.
 *                              If it is zero we probably do not detect points even if they are inside
 *                              due to rounding errors.
 */
void
t8_forest_element_points_inside (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                 const double *points, int num_points, int *is_inside, const double tolerance);

/** Find the owner process of a given element.
 * \param [in]    forest  The forest.
 * \param [in]    gtreeid The global id of the tree in which the element lies.
 * \param [in]    element The element to look for.
 * \param [in]    eclass  The element class of the tree \a gtreeid.
 * \return                The mpirank of the process that owns \a element.
 * \note The element must not exist in the forest, but an ancestor of its first
 *       descendant has to. If the element's owner is not unique, the owner of the element's
 *       first descendant is returned.
 * \note \a forest must be committed before calling this function.
 * \see t8_forest_element_find_owner_ext
 * \see t8_forest_element_owners_bounds
 */
int
t8_forest_element_find_owner (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass);

/* TODO: if set level and partition/adapt/balance all give NULL, then
 * refine uniformly and partition/adapt/balance the uniform forest. */
/** Build a uniformly refined forest on a coarse mesh.
 * \param [in]      cmesh     A coarse mesh.
 * \param [in]      scheme    An eclass scheme.
 * \param [in]      level     An initial uniform refinement level.
 * \param [in]      do_face_ghost If true, a layer of ghost elements is created for the forest.
 * \param [in]      comm      MPI communicator to use.
 * \return                    A uniform forest with coarse mesh \a cmesh, eclass_scheme
 *                            \a scheme and refinement level \a level.
 * \note This is equivalent to calling \ref t8_forest_init, \ref t8_forest_set_cmesh,
 * \ref t8_forest_set_scheme, \ref t8_forest_set_level, and \ref t8_forest_commit.
 */
t8_forest_t
t8_forest_new_uniform (t8_cmesh_t cmesh, const t8_scheme_c *scheme, const int level, const int do_face_ghost,
                       sc_MPI_Comm comm);

/** Build a adapted forest from another forest.
 * \param [in]    forest_from The forest to refine
 * \param [in]    adapt_fn    Adapt function to use
 * \param [in]    recursive   If true adaptation is recursive
 * \param [in]    do_face_ghost If true, a layer of ghost elements is created for the forest.
 * \param [in]    user_data   If not NULL, the user data pointer of the forest is set to this value.
 * \return        A new forest that is adapted from \a forest_from.
 * \note This is equivalent to calling \ref t8_forest_init, \ref t8_forest_set_adapt,
 * \ref t8_forest_set_ghost, and \ref t8_forest_commit
 */
/* TODO: make user_data const. */
t8_forest_t
t8_forest_new_adapt (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn, int recursive, int do_face_ghost,
                     void *user_data);

/** Increase the reference counter of a forest.
 * \param [in,out] forest       On input, this forest must exist with positive
 *                              reference count.  It may be in any state.
 */
void
t8_forest_ref (t8_forest_t forest);

/** Decrease the reference counter of a forest.
 * If the counter reaches zero, this forest is destroyed.
 * In this case, the forest dereferences its cmesh and scheme members.
 * \param [in,out] pforest      On input, the forest pointed to must exist
 *                              with positive reference count.  It may be in
 *                              any state.  If the reference count reaches
 *                              zero, the forest is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the forest is not modified in other ways.
 */
void
t8_forest_unref (t8_forest_t *pforest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GENERAL_H */
