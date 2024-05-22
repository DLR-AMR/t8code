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

/** \file t8_forest_private.h
 * We define routines for a forest of elements that are not part
 * of the official t8_forest.h interface but used internally.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_private */

#ifndef T8_FOREST_PRIVATE_H
#define T8_FOREST_PRIVATE_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();

/* TODO: document */

/** Check whether or not \a elements contains a (in)complete family and 
 *  return the size of it or zero if no family is considered.
 * \param [in]      forest          The forest.
 * \param [in]      ltree_id        The index of considered local tree.
 * \param [in]      el_considered   The local id of the first element in 
 *                                  \a elements in the local tree of the forest.
 * \param [in]      tscheme         The scheme for the local tree.
 * \param [in]      elements        Array of elements to consider.
 * \param [in]      elements_size   Number of elements in \a elements.
 * \return          Size of family or zero, if \a elements does not contain a
 *                  family. 
 * \note            The check works for complete and incomplete forests.
 *                  In the case of complete forests, the scheme based element 
 *                  function \see t8_element_is_family is recommended.
 * \note            If the element with id \a el_considered is not the first
 *                  family member, return 0. Therefore, if return is x > 0, 
 *                  the first x elements in \a elements form a family.
 */
int
t8_forest_is_incomplete_family (const t8_forest_t forest, const t8_locidx_t ltree_id, const t8_locidx_t el_considered,
                                t8_eclass_scheme_c *tscheme, t8_element_t **elements, const int elements_size);

/* For each tree in a forest compute its first and last descendant */
void
t8_forest_compute_desc (t8_forest_t forest);

/* Create the elements on this process given a uniform partition
 * of the coarse mesh. */
void
t8_forest_populate (t8_forest_t forest);

/** Return the eclass scheme of a given element class associated to a forest.
 * This function does not check whether the given forest is committed, use with
 * caution and only if you are sure that the eclass_scheme was set.
 * \param [in]      forest     A nearly committed forest.
 * \param [in]      eclass     An element class.
 * \return          The eclass scheme of \a eclass associated to forest.
 * \see t8_forest_set_scheme
 * \note  The forest is not required to have trees of class \a eclass.
 */
t8_eclass_scheme_c *
t8_forest_get_eclass_scheme_before_commit (t8_forest_t forest, t8_eclass_t eclass);

/** Compute the maximum possible refinement level in a forest.
 * This is the minimum over all maimum refinement level of the present element
 * classes.
 * \param [in,out] forest The forest.
 */
void
t8_forest_compute_maxlevel (t8_forest_t forest);

/** Compute the minimum possible uniform refinement level on a cmesh such
 * that no process is empty.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  scheme      The element scheme for which refinement is considered.
 * \return                  The smallest refinement level l, such that a
 *                          uniform level \a l refined forest would have no empty
 *                          processes.
 * \see t8_forest_new_uniform.
 */
int
t8_forest_min_nonempty_level (t8_cmesh_t cmesh, t8_scheme_cxx_t *scheme);

/** return nonzero if the first tree of a forest is shared with a smaller
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 * \param [in]  forest    The forest.
 * \return                True if the first tree in the forest is shared with
 *                        a smaller rank. False otherwise.
 * \note \a forest must be committed before calling this function.
 */
int
t8_forest_first_tree_shared (t8_forest_t forest);

/** return nonzero if the last tree of a forest is shared with a bigger
 * process.
 * This is the case if and only if the first descendant of the first tree that we store is
 * not the first possible descendant of that tree.
 * \param [in]  forest    The forest.
 * \return                True if the last tree in the forest is shared with
 *                        a bigger rank. False otherwise.
 * \note \a forest must be committed before calling this function.
 */
int
t8_forest_last_tree_shared (t8_forest_t forest);

/* Allocate memory for trees and set their values as in from.
 * For each tree allocate enough element memory to fit the elements of from.
 * If copy_elements is true, copy the elements of from into the element memory.
 * Do not copy the first and last desc for each tree, as this is done outside in commit
 */
void
t8_forest_copy_trees (t8_forest_t forest, t8_forest_t from, int copy_elements);

/** Given the local id of a tree in a forest, return the coarse tree of the
 * cmesh that corresponds to this tree, also return the neighbor information of
 * the tree.
 * \param [in]  forest     The forest.
 * \param [in]  ltreeid    The local id of a tree in the forest.
 * \param [out] face_neigh If not NULL a pointer to the trees face_neighbor
 *                             array is stored here on return.
 * \param [out] ttf        If not NULL a pointer to the trees tree_to_face
 *                             array is stored here on return.
 * \return                 The coarse tree that matches the forest tree with local
 *                         id \a ltreeid.
 * \see t8_cmesh_trees_get_tree_ext
 */
t8_ctree_t
t8_forest_get_coarse_tree_ext (t8_forest_t forest, t8_locidx_t ltreeid, t8_locidx_t **face_neigh, int8_t **ttf);

/** Given a forest whose trees are already filled with elements compute
 * the element offset of each local tree.
 * The element offset of a tree is the number of local elements of the forest
 * that live in all the trees with a smaller treeid.
 * \param [in,out]  forest    The forest.
 * \a forest does not need to be committed before calling this function, but all
 * elements must have been constructed.
 */
void
t8_forest_compute_elements_offset (t8_forest_t forest);

/** Return an element of a tree.
 * \param [in]  tree  The tree.
 * \param [in]  elem_in_tree The index of the element within the tree.
 * \return      Returns the element with index \a elem_in_tree of the
 *              element array of \a tree.
 */
t8_element_t *
t8_forest_get_tree_element (t8_tree_t tree, t8_locidx_t elem_in_tree);

/** Return the array of elements of a tree.
 * \param [in]  forest   The forest.
 * \param [in]  ltreeid  The local id of a local tree. Must be a valid local tree id.
 * \return      Returns the array of elements of the tree.
 * \a forest must be committed before calling this function.
 */
t8_element_array_t *
t8_forest_get_tree_element_array (t8_forest_t forest, t8_locidx_t ltreeid);

/** Find the owner process of a given element, deprecated version.
 * Use t8_forest_element_find_owner instead.
 * \param [in]     forest  The forest.
 * \param [in]     gtreeid The global id of the tree in which the element lies.
 * \param [in]     element The element to look for.
 * \param [in]     eclass  The element class of the tree \a gtreeid.
 * \param [in,out] all_owners_of_tree If not NULL, a sc_array of integers.
 *                         If the element count is zero then on output all owners
 *                         of the tree are stored.
 *                         If the element count is non-zero then it is assumed to
 *                         be filled with all owners of the tree.
 * \return                 The mpirank of the process that owns \a element.
 * \note The element must exist in the forest.
 * \note \a forest must be committed before calling this function.
 */
/* TODO: This finds the owner of the first descendant of element.
 *       We call this in owners_at_face where element is a descendant,
 *       add a flag that is true is element is a descendant, such that the
 *       first desc must not be created */
/* TODO: ext  version with parameters: lower_bound, upper_bound, is_desc/is_leaf
 *       is it really needed to construct the tree owners? Cant we just use the global
 *       offset array?
 */
int
t8_forest_element_find_owner_old (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass,
                                  sc_array_t *all_owners_of_tree);

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

/** Find the owner process of a given element, if bounds for the owner process are known.
 * \param [in]    forest  The forest.
 * \param [in]    gtreeid The global id of the tree in which the element lies.
 * \param [in]    element The element to look for.
 * \param [in]    eclass  The element class of the tree \a gtreeid.
 * \param [in]    lower_bound A known lower bound for the owner process.
 * \param [in]    upper_bound A known upper bound for the owner process.
 * \param [in]    guess   An initial guess for the owner. Must satisfy
 *                        \a lower_bound <= \a guess <= \a upper_bound
 * \return                The mpirank of the process that owns \a element.
 * \note If \a lower_bound = \a upper_bound, the function assumes that \a lower_bound
 *       is the owner process and immediately returns.
 * \note The owner p must satisfy \a lower_bound <= p <= \a upper_bound.
 * \note The element must not exist in the forest, but an ancestor of its first
 *       descendant has to. If the element's owner is not unique, the owner of the element's
 *       first descendant is returned.
 * \note \a forest must be committed before calling this function.
 * \see t8_forest_element_find_owner
 * \see t8_forest_element_owners_bounds
 */
int
t8_forest_element_find_owner_ext (t8_forest_t forest, t8_gloidx_t gtreeid, t8_element_t *element, t8_eclass_t eclass,
                                  int lower_bound, int upper_bound, int guess, int element_is_desc);

/** Perform a constant runtime check if a given rank is owner of a given element.
 * If the element is owned by more than one rank, then this check is only true
 * for the smallest.
 * \param [in]  forest      A forest.
 * \param [in]  element     An element of \a forest.
 * \param [in]  gtreeid     The global tree in which element is in.
 * \param [in]  eclass      The element class of the tree.
 * \param [in]  rank        An mpi rank.
 * \param [in]  element_is_desc This should be true, if \a element is its own first_descendant at
 *                          the maximum level. Must be false otherwise.
 * \return      True if and only if \a rank is the (first) owner process of \a element.
 */
int
t8_forest_element_check_owner (t8_forest_t forest, t8_element_t *element, t8_gloidx_t gtreeid, t8_eclass_t eclass,
                               int rank, int element_is_desc);

/** Find all owner processes that own descendant of a given element that
 * touch a given face. The element does not need to be a local element.
 * \param [in]     forest  The forest.
 * \param [in]     gtreeid The global id of the tree in which the element lies.
 * \param [in]     element The element to look for.
 * \param [in]     eclass  The element class of the tree \a gtreeid.
 * \param [in]     face    A face of \a element.
 * \param [in,out] owners  On input an array of integers. Its first and second entry
 *                         are taken as lower and upper bounds for the owner processes.
 *                         If empty, then no bounds are taken.
 *                         On output it stores
 *                         all owners of descendants of \a elem that touch \a face
 *                         in ascending order.
 */
void
t8_forest_element_owners_at_face (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                  t8_eclass_t eclass, int face, sc_array_t *owners);

/** Constant time algorithm to compute lower and upper bounds for the owner processes of a given element.
 * \param [in]     forest  The forest.
 * \param [in]     gtreeid The global id of the tree in which the element lies.
 * \param [in]     element The element to look for.
 * \param [in]     eclass  The element class of the tree \a gtreeid.
 * \param [in,out] lower   On input a known lower bound for the owner process,
 *                         on output a (better) bound.
 * \param [in,out] upper   On input a known upper bound for the owner process,
 *                         on output a (better) bound.
 *
 * \note If on input \a lower >= \a upper, then the bounds are not changed by this
 *        algorithm. We interpret \a lower = \a such that the owner is unique and equals \a lower.
 * \note \a forest must be committed before calling this function.
 * \see t8_forest_element_find_owner
 * \see t8_forest_element_owners_bounds
 */
void
t8_forest_element_owners_bounds (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                 t8_eclass_t eclass, int *lower, int *upper);

/** Constant time algorithm to compute lower and upper bounds for the owner
 * processes of the face leaves of a given element.
 * \param [in]     forest  The forest.
 * \param [in]     gtreeid The global id of the tree in which the element lies.
 * \param [in]     element The element to look for.
 * \param [in]     eclass  The element class of the tree \a gtreeid.
 * \param [in]     face    The face of \a element to consider.
 * \param [in,out] lower   On input a known lower bound for the owner process,
 *                         on output a (better) bound.
 * \param [in,out] upper   On input a known upper bound for the owner process,
 *                         on output a (better) bound.
 *
 * \note If on input \a lower >= \a upper, then the bounds are not changed by this
 *        algorithm. We interpret \a lower = \a such that the owner is unique and equals \a lower.
 * \note \a forest must be committed before calling this function.
 */
void
t8_forest_element_owners_at_face_bounds (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                         t8_eclass_t eclass, int face, int *lower, int *upper);

/** Find all owner processes that own descendant of a face neighbor of a
 *  given local element that touch the given face.
 * \param [in]     forest  The forest.
 * \param [in]     ltreeid The local id of the tree in which the element lies.
 * \param [in]     element The element, whose neighbor's face owners should be computed.
 * \param [in]     face    A face of \a element.
 * \param [in,out] owners  On input an array of integers. Its first and second entry
 *                         are taken as lower and upper bounds for the owner processes.
 *                         If empty, then no bounds are taken.
 *                         On output it stores all owners of descendants of the neighbor of
 *                         \a elem across \a face
 *                         that touch this face. If the neighbor element does not
 *                         exist, owners will be empty.
 * This is equivalent to calling t8_forest_element_face_neighbor and
 * t8_forest_element_owners_at_face for the resulting neighbor.
 * \note \a forest must be committed before calling this function.
 */
void
t8_forest_element_owners_at_neigh_face (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                        sc_array_t *owners);

/** Constant time algorithm to find bounds for the owner processes
 *  that own descendant of a face neighbor of a
 *  given local element that touch the given face.
 * \param [in]     forest  The forest.
 * \param [in]     ltreeid The local id of the tree in which the element lies.
 * \param [in]     element The element, whose neighbor's face owners should be computed.
 * \param [in]     face    A face of \a element.
 * \param [in,out] lower   On input a known lower bound for the owner process,
 *                         on output a (better) bound.
 * \param [in,out] upper   On input a known upper bound for the owner process,
 *                         on output a (better) bound.
 *
 * \note If on input \a lower >= \a upper, then the bounds are not changed by this
 *        algorithm. We interpret \a lower = \a such that the owner is unique and equals \a lower.
 * \note \a forest must be committed before calling this function.
 * This is equivalent to calling t8_forest_element_face_neighbor and
 * t8_forest_element_owners_at_face_bounds for the resulting neighbor.
 */
void
t8_forest_element_owners_at_neigh_face_bounds (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                               int face, int *lower, int *upper);

/** Construct all face neighbors of half size of a given element.
 * \param [in]     forest  The forest.
 * \param [in]     ltreeid The local tree id of the tree in which the element is.
 * \param [in]     elem    The element of which to construct the neighbors.
 * \param [in,out] neighs An array of allocated elements of the correct element class.
 *                        On output the face neighbors of \a elem across \a face of one
 *                        bigger refinement level are stored.
 * \param [in]     neigh_scheme The eclass scheme of the neighbors.
 * \param [in]     face    The number of the face of \a elem.
 * \param [in]     num_neighs The number of allocated element in \a neighs. Must match the
 *                         number of face neighbors of one bigger refinement level.
 * \param [out]    dual_face If not NULL, on output the face id's of the neighboring elements' faces.
 * \return                 The global id of the tree in which the neighbors are.
 *        -1 if there exists no neighbor across that face.
 */
t8_gloidx_t
t8_forest_element_half_face_neighbors (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *elem,
                                       t8_element_t *neighs[], t8_eclass_scheme_c *neigh_scheme, int face,
                                       int num_neighs, int dual_faces[]);

/** Iterate over all leaves of a forest and for each face compute the face neighbor
 * leaves with \ref t8_forest_leaf_face_neighbors and print their local element ids.
 * This function is meant for debugging only.
 * \param [in]    forest The forest.
 * \note Currently \a forest must be balanced.
 * \note \a forest must be committed before calling this function.
 */
void
t8_forest_print_all_leaf_neighbors (t8_forest_t forest);

/** Compute whether for a given element there exist leaf or ghost leaf elements in
 * the local forest that are a descendant of the element but not the element itself
 * \param [in]  forest    The forest.
 * \param [in]  gtreeid   The global id of the tree the element is in
 * \param [in]  element   The element
 * \param [in]  ts        The eclass scheme of \a element.
 * \return                True if in the forest there exists a local leaf or ghost
 *                        leaf that is a descendant of \a element but not equal to \a element.
 * \note If no ghost layer was created for the forest, only local elements are tested.
 * \note \a forest must be committed before calling this function.
 */
int
t8_forest_element_has_leaf_desc (t8_forest_t forest, t8_gloidx_t gtreeid, const t8_element_t *element,
                                 t8_eclass_scheme_c *ts);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_PRIVATE_H */
