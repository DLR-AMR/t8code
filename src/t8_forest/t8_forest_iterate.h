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

typedef int (*t8_forest_iterate_face_fn) (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                          int face, void *user_data, t8_locidx_t tree_leaf_index);

/*
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree
 * \param[in] element             the element for which the query is executed
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest that are descendants of \a element (or the element 
 *                                itself if \a is_leaf is true)
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 */
typedef int (*t8_forest_search_fn) (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                    const int is_leaf, const t8_element_array_t *leaf_elements,
                                    const t8_locidx_t tree_leaf_index);

/*
 * \param[in] forest              the forest
 * \param[in] ltreeid             the local tree id of the current tree
 * \param[in] element             the element for which the query is executed
 * \param[in] is_leaf             true if and only if \a element is a leaf element
 * \param[in] leaf_elements       the leaf elements in \a forest that are descendants of \a element (or the element 
 *                                itself if \a is_leaf is true)
 * \param[in] tree_leaf_index     the local index of the first leaf in \a leaf_elements
 * \param[in] queries             if not NULL, a query that is passed through from the search function
 * \param[in] query_indices       if \a query is not NULL the indices of \a query in the \a queries array from \ref 
 *                                t8_forest_search
 * \param[in, out] query_matches  if \a query is not NULL: true at the i-th index if and only if the element 'matches' 
 *                                the query of the i-th query index. 
 *                                if \a query is NULL: true if and only if the search should continue with the
 *                                children of \a element and the queries should be performed for this element.
 * \param[in] num_active_queries  The number of currently active queries. Does not have to be equal to query->elem_count,
 *                                since some queries might have been deactivated from previous calls
 */
typedef void (*t8_forest_query_fn) (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                    const int is_leaf, const t8_element_array_t *leaf_elements,
                                    const t8_locidx_t tree_leaf_index, sc_array_t *query, sc_array_t *query_indices,
                                    int *query_matches, const size_t num_active_queries);

T8_EXTERN_C_BEGIN ();

/* TODO: Document */
void
t8_forest_split_array (const t8_element_t *element, t8_element_array_t *leaf_elements, size_t *offsets);

/* TODO: comment */
/* Iterate over all leaves of an element that touch a given face of the element */
/* Callback is called in each recursive step with element as input.
 * leaf_index is only not negative if element is a leaf, in which case it indicates
 * the index of the leaf in the leaves of the tree. If it is negative, it is
 * - (index + 1) */
/* Top-down iteration and callback is called on each intermediate level.
 * If it returns false, the current element is not traversed further */
void
t8_forest_iterate_faces (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                         t8_element_array_t *leaf_elements, void *user_data, t8_locidx_t tree_lindex_of_first_leaf,
                         t8_forest_iterate_face_fn callback);

/* Perform a top-down search of the forest, executing a callback on each
 * intermediate element. The search will enter each tree at least once.
 * If the callback returns false for an element, its descendants
 * are not further searched.
 * To pass user data to the search_fn function use \ref t8_forest_set_user_data
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

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_ITERATE_H */
