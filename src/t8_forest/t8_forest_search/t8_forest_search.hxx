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

/** \file We declare a modern C++ interface for the search functionality. */

#ifndef T8_FOREST_SEARCH_HXX
#define T8_FOREST_SEARCH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <functional>

/*
 *   Discussion about C++ callback handling https://stackoverflow.com/questions/2298242/callback-functions-in-c
 *   We decided for option 4, using std::function together with templates.
*/

template <typename Udata = void>
using t8_search_element_callback
  = std::function<int (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                       const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, Udata &user_data)>;

template <typename Query_T, typename Udata = void>
using t8_search_queries_callback = std::function<void (
  t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, std::vector<Query_T> &queries,
  std::vector<int> &active_query_indices, std::vector<int> &query_matches, Udata &user_data)>;

template <typename Udata = void>
class search {
 public:
  search (t8_search_element_callback<Udata> element_callback, const t8_forest_t forest = nullptr);

  void
  update_forest (const t8_forest_t forest);
  void
  update_user_data (const Udata &udata);

  ~search ();

  void
  do_search ();

 private:
  void
  search_tree ();
  void
  search_recursion ();

  t8_search_element_callback<Udata> element_callback;

  const t8_forest_t &forest;
  const Udata &user_data;
};

template <typename Query_T, typename Udata = void>
class search_with_queries: public search<Udata> {
 public:
  search_with_queries (t8_search_element_callback<Udata> element_callback,
                       t8_search_queries_callback<Query_T, Udata> queries_callback, std::vector<Query_T> &queries,
                       const t8_forest_t forest = nullptr);

  ~search_with_queries ();

 private:
  t8_search_queries_callback<Query_T, Udata> queries_callback;
  const std::vector<Query_T> &queries;
};

#if 0
//General shape of search

template <typename Query_T>
void
do_search () (
    // init queries
    for each tree {
        search_tree {
            // init NCA and get elements
            element = NCA;
            search_recursion(element) {
                if (stop_due_to_queries) {
                    return;
                }
                ret = element_callback(element);
                if (!ret) return;
                do_queries ();
                
                // Prepare recursion
                for all children search_recursion (child[i]);
            }
            // delete NCA and elements
        }
    }
    // delete queries
}

search::do_queries ()
{
    return;
}

search::stop_due_to_queries () {
    return false;
}
 search_with_queries::do_queries () 
{
    if (queries) {
        if (!leaf) // Init new query array
        query_callback (element, queries<T>, query_matches);
        if (!leaf) // Fill new query array
    }
}
 search_with_queries::stop_due_to_queries () {
  if (queries != NULL && num_active == 0) {
    /* There are no queries left. We stop the recursion */
    return true; 
  }
}
#endif

#endif  // T8_FOREST_SEARCH_HXX
