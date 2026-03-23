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

/**
 * A call-back function used by \ref t8_forest_init_search for searching elements. Is called on an element and the
 * search criterion should be checked on that element. Return true if the search criterion is met, false otherwise.
 * 
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree in the cmesh.
 * \param[in] element             the element for which the search criterion is checked
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param[in] user_data           a user data pointer that can be set by the user
 * \returns                       non-zero if the search criterion is met, zero otherwise.
 */
typedef int (*t8_search_element_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                     const t8_element_t *element, const int is_leaf,
                                                     const t8_element_array_t *leaf_elements,
                                                     const t8_locidx_t tree_leaf_index, void *user_data);

/** A call-back function used by \ref t8_forest_init_search_with_queries for searching elements and
 * executing queries. Is called on an element and all queries are checked on that element. All positive
 * queries are passed further down to the children of the element up to leaf elements of the tree.
 * The results of the check are stored in \a query_matches.
 * 
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree in the cmesh.
 * \param[in] element             the element for which the search criterion is checked
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param[in] queries             a pointer to an array of queries
 * \param[in] user_data           a user data pointer that can be set by the user
 */
typedef int (*t8_search_queries_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                     const t8_element_t *element, const int is_leaf,
                                                     const t8_element_array_t *leaf_elements,
                                                     const t8_locidx_t tree_leaf_index, void *queries, void *user_data);

/** A call-back function used by \ref t8_forest_init_search_with_batched_queries for searching elements and
 * executing batched queries. Is called on an element and all queries are checked on that element. All positive
 * queries are passed further down to the children of the element up to leaf elements of the tree.
 * The results of the check are stored in \a query_matches.
 * 
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree in the cmesh.
 * \param[in] element             the element for which the search criterion is checked
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param[in] queries             a pointer to an array of queries
 * \param[in] active_query_indices a pointer to an array of indices of active queries in \a queries
 * \param[in,out] query_matches   a pointer to an array of length \a num_active_queries.
 *                                If query_matches[i] is true, then the element 'matches' the query of the
 *                                active query with index active_query_indices[i].
 * \param[in] user_data           a user data pointer that can be set by the user
 */
typedef void (*t8_search_batched_queries_callback_c_wrapper) (t8_forest_t forest, const t8_locidx_t ltreeid,
                                                              const t8_element_t *element, const int is_leaf,
                                                              const t8_element_array_t *leaf_elements,
                                                              const t8_locidx_t tree_leaf_index, const void *queries,
                                                              const size_t *active_query_indices, int *query_matches,
                                                              void *user_data);

/** A wrapper around the forest search context */
typedef struct t8_forest_c_search *t8_forest_search_c_wrapper;

/** Initialize the forest search context
 * \param[in,out] search           the search context to initialize
 * \param[in] element_callback     the callback function that checks the search criterion on an element
 * \param[in] forest               the forest on which the search is performed
 */
void
t8_forest_init_search (t8_forest_search_c_wrapper search, t8_search_element_callback_c_wrapper element_callback,
                       const t8_forest_t forest);

/** Update the forest on which the search is performed
 * \param[in,out] search           the search context to update
 * \param[in] forest               the new forest to use
 */
void
t8_forest_search_update_forest (t8_forest_search_c_wrapper search, const t8_forest_t forest);

/** Update the user data pointer in the search context
 * \param[in,out] search           the search context to update
 * \param[in] udata                the new user data pointer to use
 */
void
t8_forest_search_update_user_data (t8_forest_search_c_wrapper search, void *udata);

/** Perform the search
 * \param[in,out] search           the search context to use
 */
void
t8_forest_search_do_search (t8_forest_search_c_wrapper search);

/** Destroy the search context
 * \param[in,out] search           the search context to destroy
 */
void
t8_forest_search_destroy (t8_forest_search_c_wrapper search);

/**
 * A wrapper around the forest search with queries context
 */
typedef struct t8_forest_search_with_queries *t8_forest_search_with_queries_c_wrapper;

/**
 * Initialize the forest search with queries context
 * \param[in,out] search_with_queries the search with queries context to initialize
 * \param[in] element_callback        the callback function that checks the search criterion on an element
 * \param[in] queries_callback        the callback function that checks the queries on an element
 * \param[in] queries                 a pointer to an array of queries
 * \param[in] num_queries             the number of queries in the array
 * \param[in] forest                  the forest on which the search is performed
 */
void
t8_forest_init_search_with_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                    t8_search_element_callback_c_wrapper element_callback,
                                    t8_search_queries_callback_c_wrapper queries_callback, void **queries,
                                    const size_t num_queries, const t8_forest_t forest);

/** Update the forest on which the search with queries is performed
 * \param[in,out] search_with_queries the search with queries context to update
 * \param[in] forest                  the new forest to use
 */
void
t8_forest_search_with_queries_update_forest (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                             const t8_forest_t forest);

/** Update the user data pointer in the search with queries context
 * \param[in,out] search_with_queries the search with queries context to update
 * \param[in] udata                   the new user data pointer to use
 */
void
t8_forest_search_with_queries_update_user_data (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                                void *udata);

/** Update the queries in the search with queries context
 * \param[in,out] search_with_queries the search with queries context to update
 * \param[in] queries                 a pointer to an array of queries
 * \param[in] num_queries             the number of queries in the array
 */
void
t8_forest_search_with_queries_update_queries (t8_forest_search_with_queries_c_wrapper search_with_queries,
                                              void **queries, const size_t num_queries);

/** Destroy the search with queries context
 * \param[in,out] search the search with queries context to destroy
 */
void
t8_forest_search_with_queries_destroy (t8_forest_search_with_queries_c_wrapper search);

/** Perform the search with queries
 * \param[in,out] search             the search with queries context to use
 */
void
t8_forest_search_with_queries_do_search (t8_forest_search_with_queries_c_wrapper search);

/**
 * A wrapper around the forest search with batched queries context
 */
typedef struct t8_forest_search_with_batched_queries *t8_forest_search_with_batched_queries_c_wrapper;

/**
 * Initialize the forest search with batched queries context
 * \param[in,out] search_with_queries the search with batched queries context to initialize
 * \param[in] element_callback        the callback function that checks the search criterion on an element
 * \param[in] queries_callback        the callback function that checks the batched queries on an element
 * \param[in] queries                 a pointer to an array of queries
 * \param[in] num_queries             the number of queries in the array
 * \param[in] forest                  the forest on which the search is performed
 */
void
t8_forest_init_search_with_batched_queries (t8_forest_search_with_batched_queries_c_wrapper search_with_queries,
                                            t8_search_element_callback_c_wrapper element_callback,
                                            t8_search_batched_queries_callback_c_wrapper queries_callback,
                                            void **queries, const size_t num_queries, const t8_forest_t forest);

/** Update the forest on which the search with batched queries is performed
 * \param[in,out] search_with_queries the search with batched queries context to update
 * \param[in] forest                  the new forest to use
 */
void
t8_forest_search_with_batched_queries_update_forest (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, const t8_forest_t forest);

/** Update the user data pointer in the search with batched queries context
 * \param[in,out] search_with_queries the search with batched queries context to update
 * \param[in] udata                   the new user data pointer to use
 */
void
t8_forest_search_with_batched_queries_update_user_data (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void *udata);

/** Update the queries in the search with batched queries context
 * \param[in,out] search_with_queries the search with batched queries context to update
 * \param[in] queries                 a pointer to an array of queries
 * \param[in] num_queries             the number of queries in the array
 */
void
t8_forest_search_with_batched_queries_update_queries (
  t8_forest_search_with_batched_queries_c_wrapper search_with_queries, void **queries, const size_t num_queries);

/** Destroy the search with batched queries context
 * \param[in,out] search the search with batched queries context to destroy
 */
void
t8_forest_search_with_batched_queries_destroy (t8_forest_search_with_batched_queries_c_wrapper search);

/** Perform the search with batched queries
 * \param[in,out] search             the search with batched queries context to use
 */
void
t8_forest_search_with_batched_queries_do_search (t8_forest_search_with_batched_queries_c_wrapper search);

T8_EXTERN_C_END ();

#endif  // T8_FOREST_SEARCH_C_INTERFACE_H
