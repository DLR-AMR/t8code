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

/** Opaque pointer to a forest implementation. */
typedef struct t8_forest *t8_forest_t;
typedef struct t8_tree *t8_tree_t;

T8_EXTERN_C_BEGIN ();

/* TODO: There is no user_data yet */
/* TODO: if eclass is a vertex then num_outgoing/num_incoming are always
 *       1 and it is not possible to decide whether we are rfining or coarsening.
 *       Is this an issue? */
/** Callback function prototype to replace one set of elements with another.
 *
 * This is used by the adapt routine when the elements of an existing, valid
 * forest are changed.  The callback allows the user to make changes to newly
 *
 * initialized elements before the elements that they replace are destroyed.
 *
 * \param [in] forest      the forest
 * \param [in] which_tree  the local tree containing \a outgoing and \a incoming
 * \param [in] ts          the eclass scheme of the tree
 * \param [in] num_outgoing The number of outgoing elements.
 * \param [in] outgoing     The outgoing elements: after the callback, the
 *                          user_data will be destroyed. (at the current state there is no user data)
 * \param [in] num_incoming The number of incoming elements.
 * \param [in,out] incoming The incoming elements: prior to the callback,
 *                          the user_data is allocated, and the forest_init_t callback,
 *                          if it has been provided, will be called.
 *
 * If an element is being refined, num_outgoing will be 1 and num_incoming will
 * be the number of children, and vice versa if a family is being coarsened.
 */
typedef void        (*t8_forest_replace_t) (t8_forest_t forest,
                                            t8_locidx_t which_tree,
                                            t8_eclass_scheme_c * ts,
                                            int num_outgoing,
                                            t8_element_t * outgoing[],
                                            int num_incoming,
                                            t8_element_t * incoming[]);

/** Callback function prototype to decide for refining and coarsening.
 * If the \a num_elements equals the number of children then the elements
 * form a family and we decide whether this family should be coarsened
 * or only the first element should be refined.
 * Otherwise \num_elements must equal one and we consider the first entry
 * of the element array for refinement. In this case the other entries of
 * the element array are undefined.
 * \param [in] forest      the forest
 * \param [in] which_tree  the local tree containing \a elements
 * \param [in] ts          the eclass scheme of the tree
 * \param [in] num_elements the number of entries in \a elements
 * \param [in] elements    Pointers to a family or, if second entry is NULL,
 *                         pointer to one element.
 * \return greater zero if the first entry in \a elements should be refined
 *         smaller zero if the family \a elements shall be coarsened
 *         zero else.
 */
typedef int         (*t8_forest_adapt_t) (t8_forest_t forest,
                                          t8_locidx_t which_tree,
                                          t8_eclass_scheme_c * ts,
                                          int num_elements,
                                          t8_element_t * elements[]);

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
void                t8_forest_init (t8_forest_t * pforest);

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
                                          t8_scheme_cxx_t * scheme);

/** Set the initial refinement level to be used when \b forest is commited.
 * \param [in,out] forest      The forest whose level will be set.
 * \param [in]     level       The initial refinement level of \b forest, when
 *                             it is commited.
 */
void                t8_forest_set_level (t8_forest_t forest, int level);

/** Set a forest as source for copying on commiting.
 * By default, the forest takes ownership of the source \b from such that it will
 * be destroyed on calling \ref t8_forest_commit.  To keep ownership of \b
 * from, call \ref t8_forest_ref before passing it into this function.
 * This means that it is ILLEGAL to continue using \b from or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 */
void                t8_forest_set_copy (t8_forest_t forest,
                                        const t8_forest_t from);

/** Set a source forest with an adapt function to be adapted on commiting.
 * By default, the forest takes ownership of the source \b set_from such that it will
 * be destroyed on calling \ref t8_forest_commit.  To keep ownership of \b
 * set_from, call \ref t8_forest_ref before passing it into this function.
 * This means that it is ILLEGAL to continue using \b set_from or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest   The forest
 * \param [in] set_from     The source forest from which \b forest will be adapted.
 *                          We take ownership. This can be prevented by
 *                          referencing \b set_from.
 * \param [in] adapt_fn     The adapt function used on commiting.
 * \param [in] replace_fn   The replace function to be used in \b adapt_fn.
 * \param [in] recursive    A flag specifying whether adaptation is to be done recursively6
 *                          or not. If the value is zero, adaptation is not recursive
 *                          and it is recursive otherwise.
 */
void                t8_forest_set_adapt (t8_forest_t forest,
                                         const t8_forest_t set_from,
                                         t8_forest_adapt_t adapt_fn,
                                         t8_forest_replace_t replace_fn,
                                         int recursive);

/** Set the user data of a forest. This can i.e. be used to pass user defined
 * arguments to the adapt routine.
 * \param [in,out] forest   The forest
 * \param [in]     data     A pointer to user data. t8code will never touch the data.
 * The forest must not be committed before calling this function.
 * \see t8_forest_get_user_data
 */
void                t8_forest_set_user_data (t8_forest_t forest, void *data);

/** Return the user data pointer associated with a forest.
 * \param [in]     forest   The forest.
 * \return                  The user data pointer of \a forest.
 * \see t8_forest_set_user_data
 */
void               *t8_forest_get_user_data (t8_forest_t forest);

/* TODO: define weight callback function */
void                t8_forest_set_partition (t8_forest_t forest,
                                             const t8_forest_t from,
                                             int set_for_coarsening);

void                t8_forest_set_balance (t8_forest_t forest,
                                           int do_balance);
void                t8_forest_set_ghost (t8_forest_t forest, int do_ghost);

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

t8_locidx_t         t8_forest_get_num_element (t8_forest_t forest);

t8_gloidx_t         t8_forest_get_global_num_elements (t8_forest_t forest);

/** Return the element class of a forest local tree.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree in \a forest.
 * \return  The element class of the tree \a ltreeid.
 * \a forest must be committed before calling this function.
 */
t8_eclass_t         t8_forest_get_eclass (t8_forest_t forest,
                                          t8_locidx_t ltreeid);

/** Given the local id of a tree in a forest, compute the tree's local id
 * in the associated cmesh.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree in the forest.
 * \return  The local id of the tree in the cmesh associated with the forest.
 * \a forest must be committed before calling this function.
 * \note This is the inverse function of \ref t8_forest_cmesh_ltreeid_to_ltreeid.
 */
t8_locidx_t         t8_forest_ltreeid_to_cmesh_ltreeid (t8_forest_t forest,
                                                        t8_locidx_t ltreeid);

/** Given the local id of a tree in the coarse mesh of a forest, compute
 * the tree's local id in the forest.
 *  \param [in] forest    The forest.
 *  \param [in] ltreeid   The local id of a tree in the coarse mesh of \a forest.
 * \return  The local id of the tree in the forest.
 * \a forest must be committed before calling this function.
 * \note This is the inverse function of \ref t8_forest_ltreeid_to_cmesh_ltreeid.
 */
t8_locidx_t         t8_forest_cmesh_ltreeid_to_ltreeid (t8_forest_t forest,
                                                        t8_locidx_t ltreeid);

/** Given the local id of a tree in a forest, return the coarse tree of the
 * cmesh that corresponds to this tree.
 * \param [in] forest     The forest.
 * \param [in] ltreeid    The local id of a tree in the forest.
 * \return                The coarse tree that matches the forest tree with local
 *                        id \a ltreeid.
 */
t8_ctree_t          t8_forest_get_coarse_tree (t8_forest_t forest,
                                               t8_locidx_t ltreeid);

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

/** Return the number of global trees of a given forest.
 * \param [in]      forest      The forest.
 * \return          The number of global trees of that forest.
 */
t8_gloidx_t         t8_forest_get_num_global_trees (t8_forest_t forest);

/** Return a pointer to a tree in a forest.
 * \param [in]      forest      The forest.
 * \param [in]      ltree_id    The local id of the tree.
 * \return                      A pointer to the tree with local id \a ltree_id.
 * \a forest must be committed before calling this function.
 */
t8_tree_t           t8_forest_get_tree (t8_forest_t forest,
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
 */
t8_element_t       *t8_forest_get_element (t8_forest_t forest,
                                           t8_locidx_t lelement_id,
                                           t8_locidx_t * ltreeid);

/** Return the number of elements of a tree.
 * \param [in]      tree       A tree in a forest.
 * \return                     The number of elements of that tree.
 */
t8_locidx_t         t8_forest_get_tree_element_count (t8_tree_t tree);

/** Return the eclass of a tree in a forest.
 * \param [in]      forest    The forest.
 * \param [in]      ltreeid   The local id of a tree in \a forest.
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
                                                       const t8_element_t *
                                                       elem, int face);

/** Construct the face neighbor of an element, possibly across tree boundaries.
 * Returns the tree-id of the tree in which the neighbor element lies in.
 *
 * \param [in] elem The element to be considered.
 * \param [in,out] neigh On input an allocated element of the scheme of the
 *                  face_neighbors eclass.
 *                  On output, this element's data is filled with the
 *                  data of the face neighbor. If the neighbor does not exist
 *                  the data could be modified arbitrarily.
 * \param [in] face The number of the face along which the neighbor should be
 *                  constructed.
 * \return The local tree-id of the tree in which \a neigh is in.
 *        -1 if there exists no neighbor across that face.
 */
t8_locidx_t         t8_forest_element_face_neighbor (t8_forest_t forest,
                                                     t8_locidx_t ltreeid,
                                                     const t8_element_t *
                                                     elem,
                                                     t8_element_t * neigh,
                                                     int face);

/* TODO: implement */
void                t8_forest_save (t8_forest_t forest);

/** Write the forest in a parallel vtu format. There is one master
 * .pvtu file and each process writes in its own .vtu file.
 * \param [in]      forest    The forest to write.
 * \param [in]      filename  The prefix of the files where the vtk will
 *                            be stored. The master file is then filename.pvtu
 *                            and the process with rank r writes in the file
 *                            filename_r.vtu.
 * With this function the level, mpirank, treeid and element_id of each element
 * are written. For better control of the output see \ref t8_forest_vtk.h.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 */
void                t8_forest_write_vtk (t8_forest_t forest,
                                         const char *filename);

/* TODO: implement */
void                t8_forest_iterate (t8_forest_t forest);

/** Compute the coordinates of a given vertex of an element if the
 * vertex coordinates of the surrounding tree are known.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      vertices   An array storing the vertex coordinates of the tree.
 *                             It has 3*n entries, with n being the number of vertices of the tree.
 * \param [in]      corner_number The corner number of the vertex which should be computed.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the vertex.
 */
void                t8_forest_element_coordinate (t8_forest_t forest,
                                                  t8_locidx_t ltree_id,
                                                  t8_element_t * element,
                                                  const double *vertices,
                                                  int corner_number,
                                                  double *coordinates);

/** Build a uniformly refined forest on a coarse mesh.
 * \param [in]      cmesh     A coarse mesh.
 * \param [in]      scheme    An eclass scheme.
 * \param [in]      level     An initial uniform refinement level.
 * \param [in]      comm      MPI communicator to use.
 * \return                    A uniform forest with coarse mesh \a cmesh, eclass_scheme
 *                            \a scheme and refinement level \a level.
 * \note This is equivalent to calling \ref t8_forest_init, \ref t8_forest_set_cmesh,
 * \ref t8_forest_set_scheme, \ref t8_forest_set_level, and t8_forest_commit.
 */
t8_forest_t         t8_forest_new_uniform (t8_cmesh_t cmesh,
                                           t8_scheme_cxx_t * scheme,
                                           int level, sc_MPI_Comm comm);

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
void                t8_forest_unref (t8_forest_t * pforest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_H */
