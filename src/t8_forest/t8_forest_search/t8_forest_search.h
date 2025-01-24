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

/**  \file t8_forest_search.h
 * This is the C interface for the search functionality. The user can define search and query callbacks
 * to perform the search and query operations on the forest.
*/

#ifndef T8_FOREST_SEARCH_C_INTERFACE_H
#define T8_FOREST_SEARCH_C_INTERFACE_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

typedef int (*t8_search_element_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                     const t8_element_t *element, const int is_leaf,
                                                     const t8_element_array_t *leaf_elements,
                                                     const t8_locidx_t tree_leaf_index, void *user_data);

typedef int (*t8_search_queries_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                     const t8_element_t *element, const int is_leaf,
                                                     const t8_element_array_t *leaf_elements,
                                                     const t8_locidx_t tree_leaf_index, void *queries, void *user_data);

typedef void (*t8_search_batched_queries_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                              const t8_element_t *element, const int is_leaf,
                                                              const t8_element_array_t *leaf_elements,
                                                              const t8_locidx_t tree_leaf_index, const void *queries,
                                                              const size_t *active_query_indices, int *query_matches,
                                                              void *user_data);

typedef struct t8_forest_c_search *t8_forest_search_c_wrapper;

void
t8_forest_init_search (t8_forest_search_c_wrapper search, t8_search_element_callback_c_wrapper element_callback,
                       const t8_forest_t forest);
void
t8_forest_search_update_forest (t8_forest_search_c_wrapper search, const t8_forest_t forest);
void
t8_forest_search_update_user_data (t8_forest_search_c_wrapper search, void *udata);
void
t8_forest_search_do_search (t8_forest_search_c_wrapper search);
void
t8_forest_search_destroy (t8_forest_search_c_wrapper search);

typedef struct t8_forest_search_with_queries *t8_forest_search_with_queries_c_wrapper;

void
t8_forest_init_search_with_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                    t8_search_element_callback_c_wrapper element_callback,
                                    t8_search_queries_callback_c_wrapper queries_callback, void **queries,
                                    const size_t num_queries, const t8_forest_t forest);
void
t8_forest_search_with_queries_update_forest (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                             const t8_forest_t forest);
void
t8_forest_search_with_queries_update_user_data (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                                void *udata);
void
t8_forest_search_with_queries_update_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                              void **queries, const size_t num_queries);
void
t8_forest_search_with_queries_destroy (t8_forest_search_with_queries_c_wrapper search);
void
t8_forest_search_with_queries_do_search (t8_forest_search_with_queries_c_wrapper search);

typedef struct t8_forest_search_with_batched_queries *t8_forest_search_with_batched_queries_c_wrapper;

void
t8_forest_init_search_with_batched_queries (t8_forest_search_with_batched_queries_c_wrapper search_with_queries,
                                            t8_search_element_callback_c_wrapper element_callback,
                                            t8_search_batched_queries_callback_c_wrapper queries_callback,
                                            void **queries, const size_t num_queries, const t8_forest_t forest);
void
t8_forest_search_with_batched_queries_update_forest (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, const t8_forest_t forest);
void
t8_forest_search_with_batched_queries_update_user_data (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void *udata);
void
t8_forest_search_with_batched_queries_update_queries (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void **queries, const size_t num_queries);
void
t8_forest_search_with_batched_queries_destroy (t8_forest_search_with_batched_queries_c_wrapper search);
void
t8_forest_search_with_batched_queries_do_search (t8_forest_search_with_batched_queries_c_wrapper search);

T8_EXTERN_C_END ();

#endif  // T8_FOREST_SEARCH_C_INTERFACE_H
