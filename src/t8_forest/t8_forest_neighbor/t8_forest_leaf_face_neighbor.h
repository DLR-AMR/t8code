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

/** \file t8_forest_leaf_face_neighbor.h
 * Functions to query for the leaf face neighbors of a forest leaf element or leaf ghost.
 * \ref t8_forest_leaf_face_neighbors
 * \ref t8_forest_leaf_face_neighbors_ext
 * \ref t8_forest_same_level_leaf_face_neighbor_index
 */

#ifndef T8_FOREST_LEAF_FACE_NEIGHBOR_H
#define T8_FOREST_LEAF_FACE_NEIGHBOR_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

T8_EXTERN_C_BEGIN ();

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

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_LEAF_FACE_NEIGHBOR_H */
