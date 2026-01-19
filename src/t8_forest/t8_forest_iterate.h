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

/** \file t8_forest_iterate.h
 * We define various routines to iterate through (parts of) a forest and execute
 * callback functions on the leaf elements.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_iterate */

#ifndef T8_FOREST_ITERATE_H
#define T8_FOREST_ITERATE_H

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

/**
 * Callback function used in \see t8_forest_iterate_faces.
 * 
 * \param[in] forest          The forest.
 * \param[in] ltreeid         Local index of the tree containing the \a element.
 * \param[in] element         The considered element.
 * \param[in] face            The integer index of the considered face of \a element.
 * \param[in] user_data       Some user-defined data, as void pointer.
 * \param[in] tree_leaf_index Tree-local index the first leaf.
 * 
 * \return Nonzero if the element may touch the face and the top-down search shall be continued, zero otherwise.
 */
typedef int (*t8_forest_iterate_face_fn) (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                          const t8_element_t *element, const int face, void *user_data,
                                          const t8_locidx_t tree_leaf_index);

/**
 * A call-back function used by \ref t8_forest_search describing a search-criterion. Is called on an element and the 
 * search criterion should be checked on that element. Return true if the search criterion is met, false otherwise.  
 *
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree
 * \param[in] element             the element for which the search criterion is checked.
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest that are descendants of \a element (or the element 
 *                                itself if \a is_leaf is true)
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \returns                       non-zero if the search criterion is met, zero otherwise. 
 */
typedef int (*t8_forest_search_fn) (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                    const int is_leaf, const t8_element_array_t *leaf_elements,
                                    const t8_locidx_t tree_leaf_index);

/**
 * A call-back function used by \ref t8_forest_search for queries. Is called on an element and all queries are checked
 * on that element. All positive queries are passed further down to the children of the element up to leaf elements of
 * the tree. The results of the check are stored in \a query_matches. 
 * 
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree
 * \param[in] element             the element for which the queries are executed
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest that are descendants of \a element (or the element 
 *                                itself if \a is_leaf is true)
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param[in] queries             An array of queries that are checked by the function
 * \param[in] query_indices       An array of size_t entries, where each entry is an index of a query in \a queries.
 * \param[in, out] query_matches  An array of length \a num_active_queries. 
 *                                If the element is not a leave must be set to true or false at the i-th index for 
 *                                each query, specifying whether the element 'matches' the query of the i-th query 
 *                                index or not. When the element is a leaf we can return before all entries are set. 
 * \param[in] num_active_queries  The number of currently active queries (equals the number of entries of 
 *                                \a query_matches and entries of \a query_indices). 
 */
typedef void (*t8_forest_query_fn) (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                    const int is_leaf, const t8_element_array_t *leaf_elements,
                                    const t8_locidx_t tree_leaf_index, sc_array_t *queries, sc_array_t *query_indices,
                                    int *query_matches, const size_t num_active_queries);
/**
 * A call-back function used by \ref t8_forest_search_partition describing a search-criterion. Is called on an element
 * and the search criterion should be checked on that element. Return true if the search criterion is met, false
 * otherwise.
 *
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree in the cmesh. Since the cmesh has to be
 *                                replicated, it coincides with the global tree id.
 * \param[in] element             the element for which the search criterion is checked
 * \param[in] pfirst              the first processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] plast               the last processor that owns part of \a element. Guaranteed to be non-empty.
 * \returns                       non-zero if the search criterion is met, zero otherwise.
 */
typedef int (*t8_forest_partition_search_fn) (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                              const t8_element_t *element, const int pfirst, const int plast);

/**
 * A call-back function used by \ref t8_forest_search_partition for queries. Is called on an element and all queries are
 * checked on that element. All positive queries are passed further down to the children of the element. The results of
 * the check are stored in \a query_matches.
 *
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree in the cmesh. Since the cmesh has to be
 *                                replicated, it coincides with the global tree id.
 * \param[in] element             the element for which the query is executed
 * \param[in] pfirst              the first processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] plast               the last processor that owns part of \a element. Guaranteed to be non-empty.
 *                                if this is equal to \a pfirst, then the recursion will stop for
 *                                \a element's branch after this function returns.
 * \param[in] queries             an array of queries that are checked by the function
 * \param[in] query_indices       an array of size_t entries, where each entry is an index of a query in \a queries.
 * \param[in, out] query_matches  an array of length \a num_active_queries.
 *                                If the element is not a leaf must be set to true or false at the i-th index for
 *                                each query, specifying whether the element 'matches' the query of the i-th query
 *                                index or not. When the element is a leaf we can return before all entries are set.
 * \param[in] num_active_queries  The number of currently active queries (equals the number of entries of
 *                                \a query_matches and entries of \a query_indices).
 */
typedef void (*t8_forest_partition_query_fn) (const t8_forest_t forest, const t8_locidx_t ltreeid,
                                              const t8_element_t *element, const int pfirst, const int plast,
                                              void *queries, sc_array_t *query_indices, int *query_matches,
                                              const size_t num_active_queries);

T8_EXTERN_C_BEGIN ();

/** Split an array of elements according to the children of a given element E.
 *  In other words for each child C of E, find
 *  the index i, j, such that all descendants of C are
 *  elements[i], ..., elements[j-1]. 
 * 
 * \param [in] element       An element.
 * \param [in] leaf_elements An array of leaf elements of \a element. Thus, all
 *                           elements must be descendants. Sorted by linear index.
 * \param [in,out] offsets   On input an allocated array of \a num_children_of_E + 1 entries.  
 *                           On output entry i indicates the position in \a leaf_elements
 *                           where the descandents of the i-th child of E start.
 */

void
t8_forest_split_array (const t8_element_t *element, const t8_element_array_t *leaf_elements, size_t *offsets);

/**
 * Iterate over all leaves of an element that touch a given face of the element.
 * Callback is called in each recursive step with element as input.
 * leaf_index is only not negative if element is a leaf, in which case it indicates
 * the index of the leaf in the leaves of the tree. If it is negative, it is
 * - (index + 1)
 * Top-down iteration and callback is called on each intermediate level.
 * If it returns false, the current element is not traversed further 
 * 
 * \param[in] forest                    A committed forest.
 * \param[in] ltreeid                   Local index of the tree containing the \a element.
 * \param[in] element                   The considered element.
 * \param[in] face                      The integer index of the considered face of \a element.
 * \param[in] leaf_elements             The array of leaf elements that are descendants of \a element. Sorted by linear index.
 * \param[in] user_data                 The user data passed to the \a callback function.
 * \param[in] tree_lindex_of_first_leaf Tree-local index of the first leaf.
 * \param[in] callback                  The callback function.
 */
void
t8_forest_iterate_faces (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                         const int face, const t8_element_array_t *leaf_elements, void *user_data,
                         const t8_locidx_t tree_lindex_of_first_leaf, const t8_forest_iterate_face_fn callback);

/** 
 * Perform a top-down search of the forest, executing a callback on each
 * intermediate element. The search will enter each tree at least once.
 * If the callback returns false for an element, its descendants
 * are not further searched.
 * To pass user data to the search_fn function use \ref t8_forest_set_user_data.
 *  
 * \param[in] forest    The forest.
 * \param[in] search_fn The callback function describing the search criterion.
 * \param[in] query_fn  The query function.
 * \param[in] queries   The array of queries.
 */
void
t8_forest_search (t8_forest_t forest, t8_forest_search_fn search_fn, t8_forest_query_fn query_fn, sc_array_t *queries);

/** Given two forest where the elements in one forest are either direct children or
 * parents of the elements in the other forest
 * compare the two forests and for each refined element or coarsened
 * family in the old one, call a callback function providing the local indices
 * of the old and new elements.
 * \param [in]  forest_new  A forest, each element is a parent or child of an element in \a forest_old.
 * \param [in]  forest_old  The initial forest.
 * \param [in]  replace_fn  A replace callback function.
 * \note To pass a user pointer to \a replace_fn use \ref t8_forest_set_user_data
 * and \ref t8_forest_get_user_data.
 */
void
t8_forest_iterate_replace (t8_forest_t forest_new, t8_forest_t forest_old, t8_forest_replace_t replace_fn);

/**
 * Perform a top-down search of the global partition, executing a callback on
 * each intermediate element. The search will enter each tree at least once.
 * The recursion will only go down branches that are split between multiple processors.
 * This is not a collective function. It does not communicate.
 * The function expects the coarse mesh to be replicated.
 * If the callback returns false for an element, its descendants
 * are not further searched.
 * To pass user data to \b search_fn function use \ref t8_forest_set_user_data
 *
 * \param[in] forest              the forest to be searched
 * \param[in] search_fn           a search callback function called on elements
 * \param[in] query_fn            a query callback function called for all active queries of an element
 * \param[in,out] queries         an array of queries that are checked by the function
 */
void
t8_forest_search_partition (const t8_forest_t forest, t8_forest_partition_search_fn search_fn,
                            t8_forest_partition_query_fn query_fn, sc_array_t *queries);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_ITERATE_H */
