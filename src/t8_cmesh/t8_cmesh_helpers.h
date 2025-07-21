/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023, 2024 the developers

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

/** \file t8_cmesh_helpers.h
 *
 * Collection of helper routines for building cmeshes.
 */

#ifndef T8_CMESH_HELPERS_H
#define T8_CMESH_HELPERS_H

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_eclass.h>

T8_EXTERN_C_BEGIN ();

/** Sets the face connectivity information of an un-committed cmesh based on a list of tree vertices.
 * \param[in,out]   cmesh               Pointer to a t8code cmesh object. If set to NULL this argument is ignored.
 * \param[in]       ntrees              Number of coarse mesh elements resp. trees.
 * \param[in]       vertices            List of per element vertices with dimensions
 *                                      [ntrees,T8_ECLASS_MAX_CORNERS,T8_ECLASS_MAX_DIM].
 * \param[in]       eclasses            List of element classes of length [ntrees].
 * \param[in,out]   connectivity        If connectivity is not NULL the variable is filled with a pointer to an
 *                                      allocated face connectivity array. The ownership of this
 *                                      array goes to the caller. This argument is mainly used for debugging and
 *                                      testing purposes. The dimension of \a connectivity are
 *                                      [ntrees,T8_ECLASS_MAX_FACES,3].
 *                                      For each element and each face the following is stored:
 *                                      neighbor_tree_id, neighbor_dual_face_id, orientation
 * \param[in]       do_both_directions  Compute the connectivity from both neighboring sides.
 *                                      Takes much longer to compute.
 *
 * \warning  This routine might be too expensive for very large meshes. In this case, 
 *           consider to use a fully featured mesh generator.
 *
 * \note This routine does not detect periodic boundaries.
 */
void
t8_cmesh_set_join_by_vertices (t8_cmesh_t cmesh, const t8_gloidx_t ntrees, const t8_eclass_t *eclasses,
                               const double *vertices, int **connectivity, const int do_both_directions);

/** Sets the face connectivity information of an un-committed cmesh based on the cmesh stash.
 * \param[in,out]   cmesh               An uncommitted cmesh. The trees eclasses and vertices do need to be set.
 * \param[in,out]   connectivity        If connectivity is not NULL the variable is filled with a pointer to an
 *                                      allocated face connectivity array. The ownership of this
 *                                      array goes to the caller. This argument is mainly used for debugging and
 *                                      testing purposes. The dimension of \a connectivity are 
 *                                      [ntrees,T8_ECLASS_MAX_FACES,3].
 *                                      For each element and each face the following is stored:
 *                                      neighbor_tree_id, neighbor_dual_face_id, orientation
 * \param[in]       do_both_directions  Compute the connectivity from both neighboring sides. Takes much longer to compute.
 *
 * \warning This routine might be too expensive for very large meshes. In this case,
 *          consider to use a fully featured mesh generator.
 *
 * \note This routine does not detect periodic boundaries.
 */
void
t8_cmesh_set_join_by_stash (t8_cmesh_t cmesh, int **connectivity, const int do_both_directions);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_HELPERS_H */
