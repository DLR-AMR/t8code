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

#ifndef T8_FOREST_PFC_HELPER_H
#define T8_FOREST_PFC_HELPER_H
#include <t8.h>
#include <t8_data/t8_shmem.h>
#include <t8_eclass.h>

/** Make available to test */
void
t8_forest_pfc_determine_communication_range (t8_shmem_array_t partition_shmem, int rank, t8_gloidx_t relevant_begin,
                                             t8_locidx_t relevant_end, int &begin, int &end, int &num_interactions,
                                             bool skip_self);

void
t8_forest_pfc_determine_send_range (t8_shmem_array_t partition_shmem, int rank, t8_procidx_t &begin_procid,
                                    t8_procidx_t &end_procid, t8_procidx_t &num_sends, bool skip_self);

void
t8_forest_pfc_determine_recv_range (t8_shmem_array_t partition_shmem, int rank, t8_procidx_t &begin_procid,
                                    t8_procidx_t &end_procid, t8_procidx_t &num_recvs, bool skip_self);

t8_locidx_t
t8_forest_pfc_extreme_local_sibling (t8_eclass_scheme_c *scheme, t8_tree_t tree, t8_locidx_t start_element_id_in_tree,
                                     bool min_instead_max);

void
t8_forest_pfc_helper_index_in_tree_from_globalid (t8_forest_t forest, t8_gloidx_t gelement_id, t8_gloidx_t &gtree_id,
                                                  t8_eclass_scheme_c *&scheme, t8_tree_t &tree,
                                                  t8_locidx_t &index_in_tree, t8_element_t *&element);
#endif
