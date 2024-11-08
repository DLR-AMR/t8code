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
#include <numeric>

/*
 *   Discussion about C++ callback handling https://stackoverflow.com/questions/2298242/callback-functions-in-c
 *   We decided for option 4, using std::function together with templates.
*/

/**
 * \typedef t8_search_element_callback
 * \brief A callback function type used for searching elements in a forest.
 *
 * This callback function is invoked during the search process in a forest. It allows
 * custom operations to be performed on each element encountered during the search.
 *
 * \tparam Udata The type of user data passed to the callback. Defaults to void.
 *
 * \param[in] forest The forest in which the search is being performed.
 * \param[in] ltreeid The local tree ID of the current element.
 * \param[in] element A pointer to the current element being processed.
 * \param[in] is_leaf A bool indicating whether the current element is a leaf (non-zero) or not (zero).
 * \param[in] leaf_elements A pointer to an array of leaf elements.
 * \param[in] tree_leaf_index The index of the current leaf element within the tree.
 * \param[in] user_data A reference to user-defined data passed to the callback.
 *
 * \return True if the search should continue, false otherwise.
 */
template <typename Udata = void>
using t8_search_element_callback
  = std::function<bool (t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
                        const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, Udata &user_data)>;

/**
 * \typedef t8_search_queries_callback
 * \brief A callback function type used for search queries within a forest.
 *
 * \tparam Query_T The type of the query.
 * \tparam Udata The type of user data, defaults to void.
 *
 * \param[in] forest The forest in which the search is performed.
 * \param[in] ltreeid The local tree ID within the forest.
 * \param[in] element The element being queried.
 * \param[in] is_leaf A flag indicating if the element is a leaf.
 * \param[in] leaf_elements The array of leaf elements.
 * \param[in] tree_leaf_index The index of the leaf within the tree.
 * \param[in] queries A vector of queries to be processed.
 * \param[in, out] active_query_indices A vector of indices of active queries.
 * \param[in, out] query_matches A vector of query matches. Each entry corresponds to a query in the queries vector.
 * \param[in] user_data User-defined data passed to the callback.
 */
template <typename Query_T, typename Udata = void>
using t8_search_queries_callback = std::function<void (
  t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, std::vector<Query_T> &queries,
  std::vector<int> &active_query_indices, std::vector<bool> &query_matches, Udata &user_data)>;

/**
 * @brief A class to search for elements in a forest. A user-defined callback function is invoked for each element
 * to allow for custom search operations.
 * 
 * @tparam Udata 
 */
template <typename Udata = void>
class t8_search {
 public:
  /**
   * \brief Constructor of the search class. Sets the element callback, forest, and user data.
   * \param[in] element_callback The callback function to be invoked for each element during the search.
   * \param[in] forest The forest in which the search is performed.
   * \param[in] user_data The user-defined data to be passed to the callback.
   */
  t8_search (t8_search_element_callback<Udata> element_callback, const t8_forest_t forest = nullptr,
             const Udata &user_data = nullptr)
    : element_callback (element_callback), user_data (user_data)
  {
    t8_forest_ref (forest);
    this->forest = forest;
  };

  /**
   * \brief Updates the forest reference in the current object.
   *
   * This function updates the forest reference by dereferencing the current forest
   * and assigning the new forest to the object's forest member.
   *
   * \param forest The new forest to be assigned.
   */
  void
  update_forest (const t8_forest_t forest)
  {
    t8_forest_unref (&(this->forest));
    t8_forest_ref (forest);
    this->forest = forest;
  }

  /**
   * \brief Updates the user data with the provided data.
   *
   * This function sets the user data to the given value.
   *
   * \param[in] udata The new user data to be set.
   */
  void
  update_user_data (const Udata &udata)
  {
    this->user_data = udata;
  }

  /**
   * \brief Destructor for the search class.
   *
   * This destructor is responsible for unreferencing the forest object
   * associated with the search instance. It ensures that the resources
   * held by the forest object are only released if no further references exist.
   */
  ~t8_search ()
  {
    t8_forest_unref (&(this->forest));
  };

  /**
   * @brief Performs the search operation within the forest.
   *
   * This function executes the search algorithm to locate specific elements
   * within the forest structure. It performs a depth-first search on the forest, the criterion
   * to search for elements is defined by the element callback function.
   */
  void
  do_search ();

 private:
  /**
   * \brief Searches the tree for specific elements or conditions.
   *
   * This function performs a search operation on the tree.
   */
  void
  search_tree (const t8_locidx_t ltreeid);

  /**{

}
   * \brief Recursively searches for elements in the forest.
   *
   * This function performs a recursive search operation, used on each tree in the forest.
   */
  void
  search_recursion (const t8_locidx_t ltreeid, t8_element_t *element, const t8_eclass_scheme_c *ts,
                    t8_element_array_t *leaf_elements, const t8_locidx_t tree_lindex_of_first_leaf);

  bool
  stop_due_to_queries ()
  {
    return false;
  }

  void
  check_queries (const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
                 const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index)
  {
    return;
  }

  t8_search_element_callback<Udata> element_callback;

  const t8_forest_t &forest;
  const Udata &user_data;
};

template <typename Query_T, typename Udata = void>
class t8_search_with_queries: public t8_search<Udata> {
 public:
  /**
   * \brief Constructor for the search_with_queries class.
   *
   * This constructor initializes a search_with_queries object with the provided element callback,
   * queries callback, and a list of queries. It also optionally takes a forest object.
   *
   * \tparam Udata The type of user data.
   * \tparam Query_T The type of the query.
   * \param[in] element_callback A callback function for processing elements.
   * \param[in] queries_callback A callback function for processing queries.
   * \param[in] queries A vector containing the queries to be processed.
   * \param[in] forest An optional forest object. Defaults to nullptr.
   */
  t8_search_with_queries (t8_search_element_callback<Udata> element_callback,
                          t8_search_queries_callback<Query_T, Udata> queries_callback, std::vector<Query_T> &queries,
                          const t8_forest_t forest = nullptr)
    : search<Udata> (element_callback, forest), queries_callback (queries_callback), queries (queries)
  {
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  };

  /**
   * \brief Updates the list of queries with the provided queries.
   *
   * This function replaces the current list of queries with the new list
   * provided as the parameter.
   *
   * \param queries A vector containing the new queries to be set.
   */
  void
  update_queries (const std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  ~t8_search_with_queries ()
  {
    t8_forest_unref (this->forest);
  };

 private:
  bool
  stop_due_to_queries ()
  {
    return active_queries.empty ();
  }

  void
  check_queries (std::vector<size_t> new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index)
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!active_queries.empty ()) {
      std::vector<bool> query_matches (active_queries.size ());
      queries_callback (forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, queries, active_queries,
                        query_matches, user_data);
      if (!is_leaf) {
        std::for_each (active_queries.begin (), active_queries.end (), [&] (size_t iactive) {
          if (query_matches[iactive]) {
            new_active_queries.push_back (iactive);
          }
        });
      }
    }

    t8_search_queries_callback<Query_T, Udata> queries_callback;
    const std::vector<Query_T> &queries;
    std::vector<size_t> active_queries;
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
