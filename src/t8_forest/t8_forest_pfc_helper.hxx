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

/**
 * \file This file declares some helper functions used for the partition-for-coarsening feature.
*/
#ifndef T8_FOREST_PFC_HELPER_H
#define T8_FOREST_PFC_HELPER_H

#include <t8.h>
#include <t8_data/t8_shmem.h>
#include <t8_eclass.h>

/** Determine the sibling with the biggest difference in IDs (in the given direction).
 *
 * \param[in] scheme                    the refinement scheme
 * \param[in] tree                      the considered tree
 * \param[in] start_element_id_in_tree  the tree-internal ID of the considered element
 * \param[in] min_instead_max           boolean determining whether to search in direction
 *                                      of in- or decreasing IDs
 *
 * \return The extreme sibling ID within tree, i.e., the tree-internal ID of the sibling
 *         with the biggest difference to start_element_id_in_tree.
*/
t8_locidx_t
t8_forest_pfc_extreme_local_sibling (const t8_scheme_c *scheme, t8_tree_t tree, t8_locidx_t start_element_id_in_tree,
                                     bool min_instead_max);

/** Helper function for PFC that computes multiple indices for a given global element ID.
 *
 * \param[in]   forest        the forest
 * \param[in]   gelement_id   the global element ID
 * \param[out]  gtree_id      the global ID of the tree holding the element
 * \param[out]  tree          the tree holding the element
 * \param[out]  index_in_tree index of the element within the tree
 * \param[out]  element       pointer to the considered element
*/
void
t8_forest_pfc_helper_index_in_tree_from_globalid (t8_forest_t forest, t8_gloidx_t gelement_id, t8_gloidx_t &gtree_id,
                                                  t8_tree_t &tree, t8_locidx_t &index_in_tree, t8_element_t *&element);
#endif
