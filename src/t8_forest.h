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

/** \file t8_forest.h
 * We define the forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_H
#define T8_FOREST_H

#include <t8_cmesh.h>
#include <t8_element.h>
#include <t8_vtk.h>
#include <t8_data/t8_containers.h>

/** Opaque pointer to a forest implementation. */
typedef struct t8_forest *t8_forest_t;
typedef struct t8_tree *t8_tree_t;

/** This type controls, which neighbors count as ghost elements.
 * Currently, we support face-neighbors. Vertex and edge neighbors
 * will eventually be added. */
typedef enum
{
  T8_GHOST_NONE = 0,  /**< Do not create ghost layer. */
  T8_GHOST_FACES,     /**< Consider all face (codimension 1) neighbors. */
  T8_GHOST_EDGES,     /**< Consider all edge (codimension 2) and face neighbors. */
  T8_GHOST_VERTICES   /**< Consider all vertex (codimension 3) and edge and face neighbors. */
} t8_ghost_type_t;

/** This typedef is needed as a helper construct to 
 * properly be able to define a function that returns
 * a pointer to a void fun(void) function. \see t8_forest_get_user_function.
 */
typedef void        (*t8_generic_function_pointer) (void);

T8_EXTERN_C_BEGIN ();

/* TODO: if eclass is a vertex then num_outgoing/num_incoming are always
 *       1 and it is not possible to decide whether we are rfining or coarsening.
 *       Is this an issue? */
/* TODO: We may also take the local element index within the tree as parameter.
 *       Otherwise we have to search for the elements if we want pointers to them.
 */
/** Callback function prototype to replace one set of elements with another.
 *
 * This is used by the replace routine which can be called after adapt,
 * when the elements of an existing, valid
 * forest are changed. The callback allows the user to make changes to the elements
 * of the new forest that are either refined, coarsened or the same as elements in the
 * old forest.
 *
 * \param [in] forest_old      The forest that is adapted
 * \param [in] forest_new      The forest that is newly constructed from \a forest_old
 * \param [in] which_tree      The local tree containing \a outgoing and \a incoming
 * \param [in] ts              The eclass scheme of the tree
 * \param [in] num_outgoing    The number of outgoing elements.
 * \param [in] first_outgoing  The tree local index of the first outgoing element.
 *                             0 <= first_outgoing < which_tree->num_elements
 * \param [in] num_incoming    The number of incoming elements.
 * \param [in] first_incoming  The tree local index of the first incoming element.
 *                             0 <= first_incom < new_which_tree->num_elements
 *
 * If an element is being refined, num_outgoing will be 1 and num_incoming will
 * be the number of children, and vice versa if a family is being coarsened.
 * \see t8_forest_iterate_replace
 */
typedef void        (*t8_forest_replace_t) (t8_forest_t forest_old,
                                            t8_forest_t forest_new,
                                            t8_locidx_t which_tree,
                                            t8_eclass_scheme_c *ts,
                                            int num_outgoing,
                                            t8_locidx_t first_outgoing,
                                            int num_incoming,
                                            t8_locidx_t first_incoming);

/** Callback function prototype to decide for refining and coarsening.
 * If \a is_family equals 1, the first \a num_elements in \a elements
 * form a family and we decide whether this family should be coarsened
 * or only the first element should be refined.
 * Otherwise \a is_family must equal zero and we consider the first entry
 * of the element array for refinement. 
 * Entries of the element array beyond the first \a num_elements are undefined.
 * \param [in] forest       the forest to which the new elements belong
 * \param [in] forest_from  the forest that is adapted.
 * \param [in] which_tree   the local tree containing \a elements
 * \param [in] lelement_id  the local element id in \a forest_old in the tree of the current element
 * \param [in] ts           the eclass scheme of the tree
 * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements the number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 * \return greater zero if the first entry in \a elements should be refined,
 *         smaller zero if the family \a elements shall be coarsened,
 *         zero else.
 */
/* TODO: Do we really need the forest argument? Since the forest is not committed yet it
 *       seems dangerous to expose to the user. */
typedef int         (*t8_forest_adapt_t) (t8_forest_t forest,
                                          t8_forest_t forest_from,
                                          t8_locidx_t which_tree,
                                          t8_locidx_t lelement_id,
                                          t8_eclass_scheme_c *ts,
                                          const int is_family,
                                          const int num_elements,
                                          t8_element_t *elements[]);

  /** Create a new forest with reference count one.
 * This forest needs to be specialized with the t8_forest_set_* calls.
 * Currently it is manatory to either call the functions \ref
 * t8_forest_set_mpicomm, \ref t8_forest_set_cmesh, and \ref t8_forest_set_scheme,
 * or to call one of \ref t8_forest_set_copy, \ref t8_forest_set_adapt, or
 * \ref t8_forest_set_partition.  It is illegal to mix these calls, or to
 * call more than one of the three latter functions
 * Then it needs to be set up with \ref t8_forest_commit.
 * \param [in,out] pforest      On input, this pointer must be non-NULL.
 *                              On return, this pointer set to the new forest.
 */
void                t8_forest_init (t8_forest_t *pforest);

/** Check whether a forest is not NULL, initialized and not committed.
 * In addition, it asserts that the forest is consistent as much as possible.
 * \param [in] forest           This forest is examined.  May be NULL.
 * \return                      True if forest is not NULL,
 *                              \ref t8_forest_init has been called on it,
 *                              but not \ref t8_forest_commit.
 *                              False otherwise.
 */
int                 t8_forest_is_initialized (t8_forest_t forest);

/** Check whether a forest is not NULL, initialized and committed.
 * In addition, it asserts that the forest is consistent as much as possible.
 * \param [in] forest           This forest is examined.  May be NULL.
 * \return                      True if forest is not NULL and
 *                              \ref t8_forest_init has been called on it
 *                              as well as \ref t8_forest_commit.
 *                              False otherwise.
 */
int                 t8_forest_is_committed (t8_forest_t forest);

/** Check whether two committed forests have the same local elements.
 * \param [in] forest_a The first forest.
 * \param [in] forest_b The second forest.
 * \return              True if \a forest_a and \a forest_b do have the same
 *                      number of local trees and each local tree has the same
 *                      elements, that is \ref t8_element_compare returns false
 *                      for each pair of elements of \a forest_a and \a forest_b.
 * \note This function is not collective. It only returns the state on the current
 * rank.
 */
int                 t8_forest_is_equal (t8_forest_t forest_a,
                                        t8_forest_t forest_b);

/** Set the cmesh associated to a forest.
 * By default, the forest takes ownership of the cmesh such that it will be
 * destroyed when the forest is destroyed.  To keep ownership of the cmesh,
 * call \ref t8_cmesh_ref before passing it to \ref t8_forest_set_cmesh.
 * This means that it is ILLEGAL to continue using cmesh or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest       The forest whose cmesh variable will be set.
 * \param [in]     cmesh        The cmesh to be set.  We take ownership.
 *                              This can be prevented by referencing \b cmesh.
 */
void                t8_forest_set_cmesh (t8_forest_t forest,
                                         t8_cmesh_t cmesh, sc_MPI_Comm comm);

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
void                t8_forest_set_scheme (t8_forest_t forest,
                                          t8_scheme_cxx_t *scheme);

/** Set the initial refinement level to be used when \b forest is commited.
 * \param [in,out] forest      The forest whose level will be set.
 * \param [in]     level       The initial refinement level of \b forest, when
 *                             it is commited.
 * \note This setting cannot be combined with any of the derived forest methods
 * (\ref t8_forest_set_copy, \ref t8_forest_set_adapt, \ref t8_forest_set_partition,
 * and \ref t8_forest_set_balance) and overwrites any of these settings.
 * If this function is used, then the forest is created from scratch as a uniform
 * refinement of the specified cmesh (\ref t8_forest_set_cmesh, \ref t8_forest_set_scheme).
 */
void                t8_forest_set_level (t8_forest_t forest, int level);

/** Set a forest as source for copying on commiting.
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
void                t8_forest_set_copy (t8_forest_t forest,
                                        const t8_forest_t from);

/** Set a source forest with an adapt function to be adapted on commiting.
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
 * \param [in] adapt_fn     The adapt function used on commiting.
 * \param [in] recursive    A flag specifying whether adaptation is to be done recursively
 *                          or not. If the value is zero, adaptation is not recursive
 *                          and it is recursive otherwise.
 * \note This setting can be combined with \ref t8_forest_set_partition and \ref
 * t8_forest_set_balance. The order in which these operations are executed is always
 * 1) Adapt 2) Balance 3) Partition
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
/* TODO: make recursive flag to int specifying the number of recursions? */
void                t8_forest_set_adapt (t8_forest_t forest,
                                         const t8_forest_t set_from,
                                         t8_forest_adapt_t adapt_fn,
                                         int recursive);

/** Set the user data of a forest. This can i.e. be used to pass user defined
 * arguments to the adapt routine.
 * \param [in,out] forest   The forest
 * \param [in]     data     A pointer to user data. t8code will never touch the data.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_get_user_data
 */
void                t8_forest_set_user_data (t8_forest_t forest, void *data);

/** Return the user data pointer associated with a forest.
 * \param [in]     forest   The forest.
 * \return                  The user data pointer of \a forest.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_set_user_data
 */
void               *t8_forest_get_user_data (t8_forest_t forest);

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
void                t8_forest_set_user_function (t8_forest_t forest,
                                                 t8_generic_function_pointer
                                                 functrion);

/** Return the user function pointer associated with a forest.
 * \param [in]     forest   The forest.
 * \return                  The user function pointer of \a forest.
 * The forest does not need be committed before calling this function.
 * \see t8_forest_set_user_function
 */
t8_generic_function_pointer t8_forest_get_user_function (t8_forest_t forest);

/** Set a source forest to be partitioned during commit.
 * The partitioning is done according to the SFC and each rank is assinged
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
 * 1) Adapt 2) Balance 3) Partition
 * If \ref t8_forest_set_balance is called with the \a no_repartition parameter set as
 * false, it is not neccessary to call \ref t8_forest_set_partition additionally.
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
void                t8_forest_set_partition (t8_forest_t forest,
                                             const t8_forest_t set_from,
                                             int set_for_coarsening);

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
 *                          neccessary.
 * \note This setting can be combined with \ref t8_forest_set_adapt and \ref
 * t8_forest_set_balance. The order in which these operations are executed is always
 * 1) Adapt 2) Balance 3) Partition.
 * \note This setting may not be combined with \ref t8_forest_set_copy and overwrites
 * this setting.
 */
void                t8_forest_set_balance (t8_forest_t forest,
                                           const t8_forest_t set_from,
                                           int no_repartition);

/** Enable or disable the creation of a layer of ghost elements.
 * On default no ghosts are created.
 * \param [in]      forest    The forest.
 * \param [in]      do_ghost  If non-zero a ghost layer will be created.
 * \param [in]      ghost_type Controls which neighbors count as ghost elements,
 *                             currently only T8_GHOST_FACES is supported. This value
 *                             is ignored if \a do_ghost = 0.
 */
void                t8_forest_set_ghost (t8_forest_t forest, int do_ghost,
                                         t8_ghost_type_t ghost_type);

/** Like \ref t8_forest_set_ghost but with the additional options to change the
 * ghost algorithm. This is used for debugging and timing the algorithm.
 * An application should almost always use \ref t8_forest_set_ghost.
 * \param [in]      ghost_version If 1, the iterative ghost algorithm for balanced forests is used.
 *                                If 2, the iterativ algorithm for unbalanced forests.
 *                                If 3, the top-down search algorithm for unbalanced forests.
 * \see t8_forest_set_ghost
 */
void                t8_forest_set_ghost_ext (t8_forest_t forest, int do_ghost,
                                             t8_ghost_type_t ghost_type,
                                             int ghost_version);

/* TODO: use assertions and document that the forest_set (..., from) and
 *       set_load are mutually exclusive. */
void                t8_forest_set_load (t8_forest_t forest,
                                        const char *filename);

/** Compute the global number of elements in a forest as the sum
 *  of the local element counts.
 *  \param [in] forest    The forest.
 */
void                t8_forest_comm_global_num_elements (t8_forest_t forest);

/** After allocating and adding properties to a forest, commit the changes.
 * This call sets up the internal state of the forest.
 * \param [in,out] forest       Must be created with \ref t8_forest_init and
 *                              specialized with t8_forest_set_* calls first.
 */
void                t8_forest_commit (t8_forest_t forest);

/** Return the maximum allowed refinement level for any element in a forest.
 * \param [in]  forest    A forest.
 * \return                The maximum level of refinement that is allowed for
 *                        an element in this forest. It is guarenteed that any tree
 *                        in \a forest can be refined this many times and it is not
 *                        allowed to refine further.
 * \a forest must be committed before calling this function.
 * For forest with a single element class (non-hybrid) maxlevel is the maximum
 * refinement level of this element class, whilst for hybrid forests the maxlevel is
 * the minimum of all maxlevels of the element classes in this forest.
 */
int                 t8_forest_get_maxlevel (t8_forest_t forest);

/** Return the number of process local elements in the forest.
  * \param [in]  forest    A forest.
  * \return                The number of elements on this process in \a forest.
 * \a forest must be committed before calling this function.
  */
t8_locidx_t         t8_forest_get_local_num_elements (t8_forest_t forest);

/** Return the number of global elements in the forest.
  * \param [in]  forest    A forest.
  * \return                The number of elements (summed over all processes) in \a forest.
 * \a forest must be committed before calling this function.
  */
t8_gloidx_t         t8_forest_get_global_num_elements (t8_forest_t forest);

/** Return the number of ghost elements of a forest.
 * \param [in]      forest      The forest.
 * \return                      The number of ghost elements stored in the ghost
 *                              structure of \a forest. 0 if no ghosts were constructed.
 *                              \see t8_forest_set_ghost
 * \a forest must be committed before calling this function.
 */
t8_locidx_t         t8_forest_get_num_ghosts (t8_forest_t forest);

/** Return the element class of a forest local tree.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree in \a forest.
 * \return  The element class of the tree \a ltreeid.
 * \a forest must be committed before calling this function.
 */
t8_eclass_t         t8_forest_get_eclass (t8_forest_t forest,
                                          t8_locidx_t ltreeid);

/** Given a global tree id compute the forest local id of this tree.
 * If the tree is a local tree, then the local id is between 0 and the number
 * of local trees. If the tree is not a local tree, a negative number is returned.
 * \param [in]      forest The forest.
 * \param [in]      gtreeid The global id of a tree.
 * \return                 The tree's local id in \a forest, if it is a local tree.
 *                         A negative number if not.
 */
t8_locidx_t         t8_forest_get_local_id (t8_forest_t forest,
                                            t8_gloidx_t gtreeid);

/** Given the local id of a tree in a forest, compute the tree's local id
 * in the associated cmesh.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree or ghost in the forest.
 * \return  The local id of the tree in the cmesh associated with the forest.
 * \a forest must be committed before calling this function.
 * \note For forest local trees, this is the inverse function of \ref t8_forest_cmesh_ltreeid_to_ltreeid.
 */
t8_locidx_t         t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest,
                                                        t8_locidx_t ltreeid);

/** Given the local id of a tree in the coarse mesh of a forest, compute
 * the tree's local id in the forest.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree in the coarse mesh of \a forest.
 * \return  The local id of the tree in the forest. -1 if the tree is not forest local.
 * \a forest must be committed before calling this function.
 * \note For forest local trees, this is the inverse function of \ref t8_forest_ltreeid_to_cmesh_ltreeid.
 */
t8_locidx_t         t8_forest_cmesh_ltreeid_to_ltreeid (t8_forest_t forest,
                                                        t8_locidx_t lctreeid);

/** Given the local id of a tree in a forest, return the coarse tree of the
 * cmesh that corresponds to this tree.
 * \param [in] forest     The forest.
 * \param [in] ltreeid    The local id of a tree in the forest.
 * \return                The coarse tree that matches the forest tree with local
 *                        id \a ltreeid.
 */
t8_ctree_t          t8_forest_get_coarse_tree (t8_forest_t forest,
                                               t8_locidx_t ltreeid);

/** Compute the leaf face neighbors of a forest.
 * \param [in]    forest  The forest. Must have a valid ghost layer.
 * \param [in]    ltreeid A local tree id.
 * \param [in]    leaf    A leaf in tree \a ltreeid of \a forest.
 * \param [out]   neighbor_leafs Unallocated on input. On output the neighbor
 *                        leafs are stored here.
 * \param [in]    face    The index of the face across which the face neighbors
 *                        are searched.
 * \param [out]   dual_face On output the face id's of the neighboring elements' faces.
 * \param [out]   num_neighbors On output the number of neighbor leafs.
 * \param [out]   pelement_indices Unallocated on input. On output the element indices
 *                        of the neighbor leafs are stored here.
 *                        0, 1, ... num_local_el - 1 for local leafs and
 *                        num_local_el , ... , num_local_el + num_ghosts - 1 for ghosts.
 * \param [out]   pneigh_scheme On output the eclass scheme of the neighbor elements.
 * \param [in]    forest_is_balanced True if we know that \a forest is balanced, false
 *                        otherwise.
 * \note If there are no face neighbors, then *neighbor_leafs = NULL, num_neighbors = 0,
 * and *pelement_indices = NULL on output.
 * \note Currently \a forest must be balanced.
 * \note \a forest must be committed before calling this function.
 */
void                t8_forest_leaf_face_neighbors (t8_forest_t forest,
                                                   t8_locidx_t ltreeid,
                                                   const t8_element_t *leaf,
                                                   t8_element_t
                                                   **pneighbor_leafs[],
                                                   int face,
                                                   int *dual_faces[],
                                                   int *num_neighbors,
                                                   t8_locidx_t
                                                   **pelement_indices,
                                                   t8_eclass_scheme_c
                                                   **pneigh_scheme,
                                                   int forest_is_balanced);

/** Exchange ghost information of user defined element data.
 * \param[in] forest       The forest. Must be committed.
 * \param[in] element_data An array of length num_local_elements + num_ghosts
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
void                t8_forest_ghost_exchange_data (t8_forest_t forest,
                                                   sc_array_t *element_data);

/** Enable or disable profiling for a forest. If profiling is enabled, runtimes
 * and statistics are collected during forest_commit.
 * \param [in,out] forest        The forest to be updated.
 * \param [in]     set_profiling If true, profiling will be enabled, if false
 *                              disabled.
 *
 * Profiling is disabled by default.
 * The forest must not be committed before calling this function.
 * \see t8_forest_print_profile
 */
void                t8_forest_set_profiling (t8_forest_t forest,
                                             int set_profiling);

/** Print the collected statistics from a forest profile.
 * \param [in]    forest        The forest.
 *
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 */
void                t8_forest_print_profile (t8_forest_t forest);

/** Get the runtime of the last call to \ref t8_forest_adapt.
 * \param [in]   forest         The forest.
 * \return                      The runtime of adapt if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_adapt
 */
double              t8_forest_profile_get_adapt_time (t8_forest_t forest);

/** Get the runtime of the last call to \ref t8_forest_partition.
 * \param [in]   forest         The forest.
 * \param [out]  procs_sent     On output the number of processes that this rank
 *                              sent elements to in partition
 *                              if profiling was activated.
 * \return                      The runtime of partition if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_partition
 */
double              t8_forest_profile_get_partition_time (t8_forest_t forest,
                                                          int *procs_sent);

/** Get the runtime of the last call to \ref t8_forest_balance.
 * \param [in]   forest         The forest.
 * \param [out]  balance_rounts On output the number of rounds in balance
 *                              if profiling was activated.
 * \return                      The runtime of balance if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_balance
 */
double              t8_forest_profile_get_balance_time (t8_forest_t forest,
                                                        int *balance_rounds);

/** Get the runtime of the last call to \ref t8_forest_create_ghosts.
 * \param [in]   forest         The forest.
 * \param [out]  ghosts_sent    On output the number of ghost elements sent to other processes
 *                              if profiling was activated.
 * \return                      The runtime of ghost if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_set_ghost
 */
double              t8_forest_profile_get_ghost_time (t8_forest_t forest,
                                                      t8_locidx_t
                                                      *ghosts_sent);

/** Get the waittime of the last call to \ref t8_forest_ghost_exchange_data.
 * \param [in]   forest         The forest.
 * \return                      The time of ghost_exchange_data that was spent waiting
 *                              for other MPI processes, if profiling was activated.
 *                              0 otherwise.
 * \a forest must be committed before calling this function.
 * \see t8_forest_set_profiling
 * \see t8_forest_ghost_exchange_data
 */
double              t8_forest_profile_get_ghostexchange_waittime (t8_forest_t
                                                                  forest);

/** Print the ghost structure of a forest. Only used for debugging. */
void                t8_forest_ghost_print (t8_forest_t forest);

/** Change the cmesh associated to a forest to a partitioned cmesh that
 * is partitioned according to the tree distribution in the forest.
 * \param [in,out]   forest The forest.
 * \param [in]       comm   The MPI communicator that is used to partition
 *                          and commit the cmesh.
 * \param [in]       set_profiling If true, profiling for the new cmesh
 *                          will be enabled. \see t8_cmesh_set_profiling, \see t8_cmesh_print_profile
 *  \see t8_cmesh.h
 */
void                t8_forest_partition_cmesh (t8_forest_t forest,
                                               sc_MPI_Comm comm,
                                               int set_profiling);

/** Return the mpi communicator associated to a forest.
 * \param [in]      forest      The forest.
 * \return                      The mpi communicator of \a forest.
 * \a forest must be committed before calling this function.
 */
sc_MPI_Comm         t8_forest_get_mpicomm (t8_forest_t forest);

/** Return the global id of the first local tree of a forest.
 * \param [in]      forest      The forest.
 * \return                      The global id of the first local tree in \a forest.
 */
t8_gloidx_t         t8_forest_get_first_local_tree_id (t8_forest_t forest);

/** Return the number of local trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of local trees of that forest.
 */
t8_locidx_t         t8_forest_get_num_local_trees (t8_forest_t forest);

/** Return the number of ghost trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of ghost trees of that forest.
 */
t8_locidx_t         t8_forest_get_num_ghost_trees (t8_forest_t forest);

/** Return the number of global trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of global trees of that forest.
 */
t8_gloidx_t         t8_forest_get_num_global_trees (t8_forest_t forest);

/** Return the global id of a local tree or a ghost tree.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     An id 0 <= \a ltreeid < num_local_trees + num_ghosts
 *                              specifying a local tree or ghost tree.
 * \return          The global id corresponding to the tree with local id \a ltreeid.
 * \a forest must be committed before calling this function.
 */
t8_gloidx_t         t8_forest_global_tree_id (t8_forest_t forest,
                                              t8_locidx_t ltreeid);

/** Return a pointer to a tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The local id of the tree.
 * \return                      A pointer to the tree with local id \a ltree_id.
 * \a forest must be committed before calling this function.
 */
t8_tree_t           t8_forest_get_tree (t8_forest_t forest,
                                        t8_locidx_t ltree_id);

/** Return a pointer to the vertex coordinates of a tree.
 * \param [in]    forest        The forest.
 * \param [in]    ltreeid       The id of a local tree.
 * \return    If stored, a pointer to the vertex coordinates of \a tree.
 *            If no coordinates for this tree are found, NULL.
 */
double             *t8_forest_get_tree_vertices (t8_forest_t forest,
                                                 t8_locidx_t ltreeid);

/** Return the array of leaf elements of a local tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The local id of a local tree of \a forest.
 * \return                      An array of t8_element_t * storing all leaf elements
 *                              of this tree.
 */
t8_element_array_t *t8_forest_tree_get_leafs (t8_forest_t forest,
                                              t8_locidx_t ltree_id);

/** Return a cmesh associated to a forest.
 * \param [in]      forest      The forest.
 * \return          The cmesh associated to the forest.
 */
t8_cmesh_t          t8_forest_get_cmesh (t8_forest_t forest);

/** Return an element of the forest.
 * \param [in]      forest      The forest.
 * \param [in]      lelement_id The local id of an element in \a forest.
 * \param [out]     ltreeid     If not NULL, on output the local tree id of the tree in which the
 *                              element lies in.
 * \return          A pointer to the element. NULL if this element does not exist.
 * \note This function performs a binary search. For constant access, use \ref t8_forest_get_element_in_tree
 * \a forest must be committed before calling this function.
 */
t8_element_t       *t8_forest_get_element (t8_forest_t forest,
                                           t8_locidx_t lelement_id,
                                           t8_locidx_t *ltreeid);

/** Return an element of a local tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     An id of a local tree in the forest.
 * \param [in]      leid_in_tree The index of an element in the tree.
 * \return          A pointer to the element.
 * \note If the tree id is know, this function should be preferred over \ref t8_forest_get_element.
 * \a forest must be committed before calling this function.
 */
t8_element_t       *t8_forest_get_element_in_tree (t8_forest_t forest,
                                                   t8_locidx_t ltreeid,
                                                   t8_locidx_t leid_in_tree);

/** Return the number of elements of a tree.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     A local id of a tree.
 * \return                      The number of elements in the local tree \a ltreeid.
 */
t8_locidx_t         t8_forest_get_tree_num_elements (t8_forest_t forest,
                                                     t8_locidx_t ltreeid);

/** Return the element offset of a local tree, that is the number of elements
 * in all trees with smaller local treeid.
 * \param [in]      forest      The forest.
 * \param [in]      ltreeid     A local id of a tree.
 * \return                      The number of leaf elements on all local tree with
 *                              id < \a ltreeid.
 * \note \a forest must be committed before calling this function.
 */
t8_locidx_t         t8_forest_get_tree_element_offset (t8_forest_t forest,
                                                       t8_locidx_t ltreeid);

/** Return the number of elements of a tree.
 * \param [in]      tree       A tree in a forest.
 * \return                     The number of elements of that tree.
 */
t8_locidx_t         t8_forest_get_tree_element_count (t8_tree_t tree);

/** Return the eclass of a tree in a forest.
 * \param [in]      forest    The forest.
 * \param [in]      ltreeid   The local id of a tree (local or ghost) in \a forest.
 * \return                    The element class of the tree with local id \a ltreeid.
 */
t8_eclass_t         t8_forest_get_tree_class (t8_forest_t forest,
                                              t8_locidx_t ltreeid);

/** Compute the global index of the first local element of a forest.
 * This function is collective.
 * \param [in]     forest       A committed forest, whose first element's index is computed.
 * \return         The global index of \a forest's first local element.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 */
t8_gloidx_t         t8_forest_get_first_local_element_id (t8_forest_t forest);

/** Return the element scheme associated to a forest.
 * \param [in]      forest.     A committed forest.
 * \return          The element scheme of the forest.
 * \see t8_forest_set_scheme
 */
t8_scheme_cxx_t    *t8_forest_get_scheme (t8_forest_t forest);

/** Return the eclass scheme of a given element class associated to a forest.
 * \param [in]      forest.     A committed forest.
 * \param [in]      eclass.     An element class.
 * \return          The eclass scheme of \a eclass associated to forest.
 * \see t8_forest_set_scheme
 * \note  The forest is not required to have trees of class \a eclass.
 */
t8_eclass_scheme_c *t8_forest_get_eclass_scheme (t8_forest_t forest,
                                                 t8_eclass_t eclass);

/** Return the eclass of the tree in which a face neighbor of a given element
 * lies.
 * \param [in]      forest.     A committed forest.
 * \param [in]      ltreeid.    The local tree in which the element lies.
 * \param [in]      elem.       An element in the tree \a ltreeid.
 * \param [in]      face.       A face number of \a elem.
 * \return                      The local tree id of the tree in which the face
 *                              neighbor of \a elem across \a face lies.
 */
t8_eclass_t         t8_forest_element_neighbor_eclass (t8_forest_t forest,
                                                       t8_locidx_t ltreeid,
                                                       const t8_element_t
                                                       *elem, int face);

/** Construct the face neighbor of an element, possibly across tree boundaries.
 * Returns the global tree-id of the tree in which the neighbor element lies in.
 *
 * \param [in] elem The element to be considered.
 * \param [in,out] neigh On input an allocated element of the scheme of the
 *                  face_neighbors eclass.
 *                  On output, this element's data is filled with the
 *                  data of the face neighbor. If the neighbor does not exist
 *                  the data could be modified arbitrarily.
 * \param [in] neigh_scheme The eclass scheme of \a neigh.
 * \param [in] face The number of the face along which the neighbor should be
 *                  constructed.
 * \param [out] neigh_face The number of the face viewed from perspective of \a neigh.
 * \return The global tree-id of the tree in which \a neigh is in.
 *        -1 if there exists no neighbor across that face.
 */
t8_gloidx_t         t8_forest_element_face_neighbor (t8_forest_t forest,
                                                     t8_locidx_t ltreeid,
                                                     const t8_element_t *elem,
                                                     t8_element_t *neigh,
                                                     t8_eclass_scheme_c
                                                     *neigh_scheme, int face,
                                                     int *neigh_face);

/* TODO: implement */
void                t8_forest_save (t8_forest_t forest);

/** Write the forest in a parallel vtu format. Extended version.
 * See \ref t8_forest_write_vtk for the standard version of this function.
 * Writes one master .pvtu file and each process writes in its own .vtu file.
 * If linked and not otherwise specified, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This may change in accordance with \a write_ghosts, \a write_curved and 
 * \a do_not_use_API, because the export of ghosts is not yet available with 
 * the VTK API and the export of curved elements is not available with the
 * inbuilt function to write ASCII files. The function will for example
 * still use the VTK API to satisfy \a write_curved, even if \a do_not_use_API 
 * is set to true.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 * \param [in]      forest              The forest to write.
 * \param [in]      fileprefix          The prefix of the files where the vtk will
 *                                      be stored. The master file is then fileprefix.pvtu
 *                                      and the process with rank r writes in the file
 *                                      fileprefix_r.vtu.
 * \param [in]      write_treeid        If true, the global tree id is written for each element.
 * \param [in]      write_mpirank       If true, the mpirank is written for each element.
 * \param [in]      write_level         If true, the refinement level is written for each element.
 * \param [in]      write_element_id    If true, the global element id is written for each element.
 * \param [in]      write_ghosts        If true, each process additionally writes its ghost elements.
 *                                      For ghost element the treeid is -1.
 * \param [in]      write_curved        If true, write the elements as curved element types from vtk.
 * \param [in]      do_not_use_API      Do not use the VTK API, even if linked and available.
 * \param [in]      num_data            Number of user defined double valued data fields to write.
 * \param [in]      data                Array of t8_vtk_data_field_t of length \a num_data
 *                                      providing the user defined per element data.
 *                                      If scalar and vector fields are used, all scalar fields
 *                                      must come first in the array.
 * \return  True if successful, false if not (process local).
 * See also \ref t8_forest_write_vtk .
 */
int                 t8_forest_write_vtk_ext (t8_forest_t forest,
                                             const char *fileprefix,
                                             int write_treeid,
                                             int write_mpirank,
                                             int write_level,
                                             int write_element_id,
                                             int write_ghosts,
                                             int write_curved,
                                             int do_not_use_API,
                                             int num_data,
                                             t8_vtk_data_field_t *data);

/** Write the forest in a parallel vtu format. Writes one master
 * .pvtu file and each process writes in its own .vtu file.
 * If linked, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This function writes the forest elements, the tree id, element level, mpirank and element id as data.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 * For more options use \ref t8_forest_write_vtk_ext
 * \param [in]      forest              The forest to write.
 * \param [in]      fileprefix          The prefix of the files where the vtk will
 *                                      be stored. The master file is then fileprefix.pvtu
 *                                      and the process with rank r writes in the file
 *                                      fileprefix_r.vtu.
 * \return  True if successful, false if not (process local).
 */
int                 t8_forest_write_vtk (t8_forest_t forest,
                                         const char *fileprefix);

/* TODO: implement */
void                t8_forest_iterate (t8_forest_t forest);

/** Compute the coordinates of a given vertex of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      corner_number The corner number, in Z-order, of the vertex which should be computed.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the vertex.
 */
void                t8_forest_element_coordinate (t8_forest_t forest,
                                                  t8_locidx_t ltree_id,
                                                  const t8_element_t *element,
                                                  int corner_number,
                                                  double *coordinates);

/** Compute the coordinates of the centroid of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * The centroid is the sum of all corner vertices divided by the number of corners.
 * The centroid can be seen as the midpoint of an element and thus can for example be used
 * to compute level-set values or the distance between two elements.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the centroid.
 */
void                t8_forest_element_centroid (t8_forest_t forest,
                                                t8_locidx_t ltreeid,
                                                const t8_element_t *element,
                                                double *coordinates);

/** Compute the diameter of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       For lines the value is exact while for other element types it is only
 *                             an approximation.
 */
double              t8_forest_element_diam (t8_forest_t forest,
                                            t8_locidx_t ltreeid,
                                            const t8_element_t *element);

/** Compute the volume of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       This function assumes d-linear interpolation for the
 *                             tree vertex coordinates.
 *                             \a forest must be committed when calling this function.
 */
double              t8_forest_element_volume (t8_forest_t forest,
                                              t8_locidx_t ltreeid,
                                              const t8_element_t *element);

/** Compute the area of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \return                     The area of \a face.
 * \a forest must be committed when calling this function.
 */
double              t8_forest_element_face_area (t8_forest_t forest,
                                                 t8_locidx_t ltreeid,
                                                 const t8_element_t *element,
                                                 int face);

/** Compute the vertex coordinates of the centroid of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the centroid of \a face.
 * \a forest must be committed when calling this function.
 */
void                t8_forest_element_face_centroid (t8_forest_t forest,
                                                     t8_locidx_t ltreeid,
                                                     const t8_element_t
                                                     *element, int face,
                                                     double centroid[3]);

/** Compute the normal vector of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the normal vector of \a element at \a face.
 * \a forest must be committed when calling this function.
 */
void                t8_forest_element_face_normal (t8_forest_t forest,
                                                   t8_locidx_t ltreeid,
                                                   const t8_element_t
                                                   *element, int face,
                                                   double normal[3]);

/** Query whether a given point lies inside an element or not. For bilinearly interpolated elements.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      point      3-dimensional coordinates of the point to check
 * \param [in]      tolerance  tolerance that we allow the point to not exactly match the element.
 *                             If this value is larger we detect more points.
 *                             If it is zero we probably do not detect points even if they are inside
 *                             due to rounding errors.
 * \return          True (non-zero) if \a point lies within \a element, false otherwise.
 *                  The return value is also true if the point lies on the element boundary.
 *                  Thus, this function may return true for different leaf elements, if they
 *                  are neighbors and the point lies on the common boundary.
 */
int                 t8_forest_element_point_inside (t8_forest_t forest,
                                                    t8_locidx_t ltreeid,
                                                    const t8_element_t
                                                    *element,
                                                    const double point[3],
                                                    const double tolerance);

/* TODO: if set level and partition/adapt/balance all give NULL, then
 * refine uniformly and partition/adapt/balance the unfiform forest. */
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
t8_forest_t         t8_forest_new_uniform (t8_cmesh_t cmesh,
                                           t8_scheme_cxx_t *scheme,
                                           int level, int do_face_ghost,
                                           sc_MPI_Comm comm);

/** Build a adapted forest from another forest.
 * \param [in]    forest_from The forest to refine
 * \param [in]    adapt_fn    Adapt function to use
 * \param [in]    replace_fn  Replace function to use
 * \param [in]    recursive   If true adptation is recursive
 * \param [in]    do_face_ghost If true, a layer of ghost elements is created for the forest.
 * \param [in]    user_data   If not NULL, the user data pointer of the forest is set to this value.
 * \return        A new forest that is adapted from \a forest_from.
 * \note This is equivalent to calling \ref t8_forest_init, \ref t8_forest_set_adapt,
 * \red t8_forest_set_ghost, and \ref t8_forest_commit
 */
/* TODO: make user_data const. */
t8_forest_t         t8_forest_new_adapt (t8_forest_t forest_from,
                                         t8_forest_adapt_t adapt_fn,
                                         int recursive, int do_face_ghost,
                                         void *user_data);

/** Increase the reference counter of a forest.
 * \param [in,out] forest       On input, this forest must exist with positive
 *                              reference count.  It may be in any state.
 */
void                t8_forest_ref (t8_forest_t forest);

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
void                t8_forest_unref (t8_forest_t *pforest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_H */
