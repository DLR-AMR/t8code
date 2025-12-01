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

/**  \file t8_forest_search.hxx 
 * A C++ interface for the search functionality. The user can define search and query callbacks
 * to perform the search and query operations on the forest. Implementation details regarding the 
 * callback handling are given by https://stackoverflow.com/questions/2298242/callback-functions-in-c
 * We decided for option 4, using std::function together with templates.
*/

#ifndef T8_FOREST_SEARCH_HXX
#define T8_FOREST_SEARCH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest.h>  // Ensure t8_forest_t is defined
#include <functional>
#include <numeric>
#include <ranges>

/**
 *  t8_search_element_callback
 *  A callback function type used for searching elements in a forest.
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
using t8_search_element_callback = std::function<bool (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, Udata *user_data)>;

/**
 *  A callback function type used for search queries within a forest.
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
 * \param[in] query A single query to be processed.
 * \param[in] user_data User-defined data passed to the callback.
 */
template <typename Query_T, typename Udata = void>
using t8_search_query_callback = std::function<bool (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, const Query_T &query, Udata *user_data)>;

/**
 *  t8_search_batched_queries_callback
 *  A callback function type used for search queries within a forest. Processes a batch of queries.
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
using t8_search_batched_queries_callback = std::function<void (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, const std::vector<Query_T> &queries,
  const std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches, Udata *user_data)>;

/**
 *  t8_partition_search_element_callback
 *  A callback function type used for searching elements in the partition of a forest.
 *
 * This callback function is invoked during the partition search process in a forest. It allows
 * custom operations to be performed on each element encountered during the search.
 *
 * \tparam Udata The type of user data passed to the callback. Defaults to void.
 *
 * \param[in] forest The forest whose partition is searched.
 * \param[in] ltreeid The local tree ID of the current tree in the cmesh.
 * \param[in] element A pointer to the current element being processed.
 * \param[in] pfirst The first processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] plast The last processor that owns part of \a element. Guaranteed to be non-empty.
 *
 * \return True, if the search should continue, false otherwise.
 */
template <typename Udata = void>
using t8_partition_search_element_callback
  = std::function<bool (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                        const int pfirst, const int plast, Udata *user_data)>;

/**
 *  t8_partition_search_query_callback
 *  A callback function type used for searching queries in the partition of a forest.
 *
 * \tparam Query_T The type of the query.
 * \tparam Udata The type of user data, defaults to void.
 *
 * \param[in] forest The forest whose partition is searched.
 * \param[in] ltreeid The local tree ID within the forest.
 * \param[in] element The element being queried.
 * \param[in] pfirst The first processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] plast The last processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] query A single query to be processed.
 * \param[in] user_data User-defined data passed to the callback.
 */
template <typename Query_T, typename Udata = void>
using t8_partition_search_query_callback
  = std::function<bool (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                        const int pfirst, const int plast, const Query_T &query, Udata *user_data)>;

/**
 *  t8_partition_search_batched_queries_callback
 *  A callback function type used for searching queries in the partition of a forest. Processes a batch of queries.
 *
 * \tparam Query_T The type of the query.
 * \tparam Udata The type of user data, defaults to void.
 *
 * \param[in] forest The forest whose partition is searched.
 * \param[in] ltreeid The local tree ID within the forest.
 * \param[in] element The element being queried.
 * \param[in] pfirst The first processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] plast The last processor that owns part of \a element. Guaranteed to be non-empty.
 * \param[in] queries A vector of queries to be processed.
 * \param[in, out] active_query_indices A vector of indices of active queries.
 * \param[in, out] query_matches A vector of query matches. Each entry corresponds to a query in the queries vector.
 * \param[in] user_data User-defined data passed to the callback.
 */
template <typename Query_T, typename Udata = void>
using t8_partition_search_batched_queries_callback = std::function<void (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const int pfirst, const int plast,
  const std::vector<Query_T> &queries, const std::vector<size_t> &active_query_indices,
  std::vector<bool> &query_matches, Udata *user_data)>;

/**
 * A base class that performs a search in a forest.
 */
class t8_search_base {
 public:
  /**   Constructor for the t8_search_base class.
   *
   *
   * This constructor initializes a t8_search_base object with the given forest.
   * If the forest is not null, it increments the reference count of the forest
   * and asserts that the forest is committed.
   *
   * \param[in] forest A pointer to a t8_forest_t object. Defaults to nullptr.
   */
  t8_search_base (t8_forest_t forest = nullptr): forest (forest)
  {
    if (forest != nullptr) {
      t8_forest_ref (forest);
      T8_ASSERT (t8_forest_is_committed (forest));
    }
  }

  /**   Update the forest for the search.
   *
   * This function updates the forest for the search. If the current forest is not null,
   * it decrements the reference count of the forest. It then asserts that the new forest
   * is not null and is committed. Finally, it increments the reference count of the new forest.
   *
   * \param[in] forest A pointer to a t8_forest_t object.
   */
  void
  update_forest (t8_forest_t forest)
  {
    if (this->forest != nullptr) {
      t8_forest_unref (&(this->forest));
    }
    T8_ASSERT (forest != nullptr);
    T8_ASSERT (t8_forest_is_committed (forest));
    t8_forest_ref (forest);
    this->forest = forest;
  }

  /**   Destructor for the t8_search_base class.
   *
   * This destructor decrements the reference count of the forest if it is not null.
   */
  ~t8_search_base ()
  {
    if (this->forest != nullptr) {
      t8_forest_unref (&(this->forest));
    }
  }

  /**   Perform the search.
   *
   * This function performs the search in the forest. 
   */
  void
  do_search ();

  /**
   * The forest on which the search is performed.
   */
  t8_forest_t forest;

 private:
  /** @brief Searches a tree within the forest.
   *
   * This function performs a search operation on a tree identified by the given local tree ID.
   * It uses the \a search_recursion function to perform the search.
   *
   * \param[in] ltreeid The local tree ID of the tree to be searched.
   */
  void
  search_tree (const t8_locidx_t ltreeid);

  /**  Recursively searches the tree.
   *
   * This function performs a recursive search operation on the tree identified by the given local tree ID.
   * It uses the given \a element_callback function to process each element encountered during the search.
   * If a query_callback function is provided, it is used to process queries during the search.
   *
   * \param[in] ltreeid The local tree ID of the tree to be searched.
   * \param[in] element The element to be searched.
   * \param[in] ts The element class scheme.
   * \param[in] leaf_elements The array of leaf elements.
   * \param[in] tree_lindex_of_first_leaf The index of the first leaf in the tree.
   */
  void
  search_recursion (const t8_locidx_t ltreeid, t8_element_t *element, const t8_scheme *ts,
                    t8_element_array_t *leaf_elements, const t8_locidx_t tree_lindex_of_first_leaf);

  /**  Checks if the search should stop due to empty queries.
   * 
   */
  virtual bool
  stop_due_to_queries ()
    = 0;

  /**  Checks an element during the search.
   *
   * This function is called for each element encountered during the search.
   * It passes the arguments to the callback function provided by the user.
   *
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] is_leaf A bool indicating whether the current element is a leaf (non-zero) or not (zero).
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current leaf element within the tree.
   *
   * \return True if the search should continue, false otherwise.
   */
  virtual bool
  check_element (const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
                 const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index)
    = 0;

  /**  Checks queries during the search.
   * 
   * This function is called to check queries during the search.
   * It passes the arguments to the callback function provided by the user.
   * 
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] is_leaf A bool indicating whether the current element is a leaf (non-zero) or not (zero).
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current leaf element within the tree.
   */
  virtual void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index)
    = 0;

  /**
   *   Function the gives the user the opportunity to update the queries after 
   *         each step in the recursion.
   *
   * \param old_query_indices 
   */
  virtual void
  update_queries (std::vector<size_t> &old_query_indices)
    = 0;

  /**
   * Function gives the user the opportunity to set the queries to the initial
   *        full set before searching each tree.
   *
   */
  virtual void
  init_queries ()
    = 0;
};

/**
 * A class that performs a search in a forest.
 * \tparam Udata The type of user data to be used in the search.
 */
template <typename Udata = void>
class t8_search: public t8_search_base {
 public:
  /**
  * Constructor for the t8_search class.
  * \param[in] element_callback A callback function to be called for each element during the search.
  * \param[in] forest A pointer to the forest to be searched. Defaults to nullptr
  * \param[in] user_data A pointer to user-defined data to be passed to the callback function. Defaults to nullptr.
  */
  t8_search (t8_search_element_callback<Udata> element_callback, t8_forest_t forest = nullptr,
             Udata *user_data = nullptr)
    : t8_search_base (forest), element_callback (element_callback)
  {
    if (user_data != nullptr) {
      this->user_data = user_data;
    }
  }

  /**  Updates the user data associated with the object.
   *
   * This function sets the user data pointer to the provided Udata object.
   *
   * \param[in] udata A pointer to the Udata object to be set as the user data.
   */
  void
  update_user_data (Udata *udata)
  {
    if (udata != nullptr) {
      this->user_data = udata;
    }
  }

  /**
   * Destructor for the t8_search class.
   */
  virtual ~t8_search () = default;

  /**
   * The user data associated with the object.
   */
  Udata *user_data;

 private:
  /**
   * Checks if the search should stop due to empty queries.
   */
  bool
  stop_due_to_queries () override
  {
    return false;
  }

  /**
   * Checks an element during the search by invoking the user-defined callback function.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] is_leaf A bool indicating whether the current element is a leaf or not.
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current leaf element within the tree.
   * \return True if the search should continue, false otherwise.
   */
  bool
  check_element (const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
                 const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index) override
  {
    T8_ASSERT (t8_forest_is_committed (this->forest));
    return this->element_callback (this->forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, user_data);
  }

  /**
   * A no-op implementation of the check_queries function.
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] is_leaf A bool indicating whether the current element is a leaf or not.
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current leaf element within the tree.
   */
  void
  check_queries ([[maybe_unused]] std::vector<size_t> &new_active_queries, [[maybe_unused]] const t8_locidx_t ltreeid,
                 [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const bool is_leaf,
                 [[maybe_unused]] const t8_element_array_t *leaf_elements,
                 [[maybe_unused]] const t8_locidx_t tree_leaf_index) override
  {
    return;
  }

  /**
   * A no-op implementation of the update_queries function.
   * \param[in] old_query_indices A vector of indices of old queries.
   */
  void
  update_queries ([[maybe_unused]] std::vector<size_t> &old_query_indices) override
  {
    return;
  }

  /**
   * A no-op implementation of the init_queries function.
   */
  void
  init_queries () override
  {
    return;
  }

  /**
   * The element callback function used during the search.
   */
  t8_search_element_callback<Udata> element_callback;
};

/** 
 *  A class that performs a search in a forest with queries.
 * Uses a filter-view to filter out the active queries. It is recommended to use this version of the search
 * if the number of queries is small or if the queries do not need any further computations to be evaluated. 
 * 
 * \tparam Query_T The type of queries
 * \tparam Udata The type of the user data, defaults to void.
*/
template <typename Query_T, typename Udata = void>
class t8_search_with_queries: public t8_search<Udata> {
 public:
  /**
   * Constructor for the t8_search_with_queries class.
   * \param[in] element_callback A callback function to be called for each element during the search.
   * \param[in] queries_callback A callback function to be called for each query during the search.
   * \param[in] queries A reference to a vector of queries to be processed.
   * \param[in] forest A pointer to the forest to be searched. Defaults to nullptr
   * \param[in] user_data A pointer to user-defined data to be passed to the callback function. Defaults to nullptr.
   */
  t8_search_with_queries (t8_search_element_callback<Udata> element_callback,
                          t8_search_query_callback<Query_T, Udata> queries_callback, std::vector<Query_T> &queries,
                          const t8_forest_t forest = nullptr, Udata *user_data = nullptr)
    : t8_search<Udata> (element_callback, forest, user_data), queries_callback (queries_callback), queries (queries)
  {
  }

  /**
   * Update the queries for the search.
   * \param[in] queries A reference to a vector of queries to be processed.
   */
  void
  update_queries (std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  /**
   * Destructor for the t8_search_with_queries class.
   */
  ~t8_search_with_queries ()
  {
  }

 private:
  /**
   * Checks if the search should stop due to empty queries.
   */
  bool
  stop_due_to_queries () override
  {
    return this->active_queries.empty ();
  }

  /**
   * Checks the queries against the current element.
   * \param[in] new_active_queries A reference to a vector of new active query indices.
   * \param[in] ltreeid The ID of the local tree.
   * \param[in] element A pointer to the current element being checked.
   * \param[in] is_leaf A boolean indicating if the current element is a leaf.
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current element in the leaf array.
   */
  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements,
                 const t8_locidx_t tree_leaf_index) override
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!this->active_queries.empty ()) {
      std::copy_if (this->active_queries.begin (), this->active_queries.end (), std::back_inserter (new_active_queries),
                    [&] (size_t &query_index) {
                      return this->queries_callback (this->forest, ltreeid, element, is_leaf, leaf_elements,
                                                     tree_leaf_index, queries[query_index], this->user_data);
                    });
    }
  }

  /**
   * Updates the queries for the search.
   * \param[in] old_query_indices A reference to a vector of queries to be processed.
   */
  void
  update_queries (std::vector<size_t> &old_query_indices) override
  {
    this->active_queries = old_query_indices;
  }

  /**
   * Initializes the queries for the search.
   */
  void
  init_queries () override
  {
    this->active_queries.resize (this->queries.size ());
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  }

  /**
   * A reference to a vector of queries to be processed.
   */
  t8_search_query_callback<Query_T, Udata> queries_callback;
  /**
   * A reference to a vector of queries to be processed.
   */
  std::vector<Query_T> &queries;
  /**
   * A vector of indices of active queries.
   */
  std::vector<size_t> active_queries;
};

/**
 * A class that performs a search in a forest with batched queries.
 * 
 * All active queries are passed to the callback function, which processes them in a batch. It is recommended to 
 * use this version of the searcch if further computations have to be done to evaluate the queries. That way these
 * precomputations are not done for every call to the callback again and only have to be evaluated once per call.
 * 
 * \tparam Query_T The type of queries
 * \tparam Udata The type of the user data, defaults to void.
 */
template <typename Query_T, typename Udata = void>
class t8_search_with_batched_queries: public t8_search<Udata> {
 public:
  /**
  * Constructor for the t8_search_with_batched_queries class.
  * \param[in] element_callback A callback function to be called for each element during the search.
  * \param[in] queries_callback A callback function to be called for processing batched queries.
  * \param[in] queries A reference to a vector of queries to be processed.
  * \param[in] forest A pointer to the forest to be searched. Defaults to nullptr
  * \param[in] user_data A pointer to user-defined data to be passed to the callback function. Defaults to nullptr.
  */
  t8_search_with_batched_queries (t8_search_element_callback<Udata> element_callback,
                                  t8_search_batched_queries_callback<Query_T, Udata> queries_callback,
                                  std::vector<Query_T> &queries, const t8_forest_t forest = nullptr,
                                  Udata *user_data = nullptr)
    : t8_search<Udata> (element_callback, forest, user_data), queries_callback (queries_callback), queries (queries)
  {
  }

  /**
   * Update the queries for the search.
   * \param[in] queries A reference to a vector of queries to be processed.
   */
  void
  update_queries (std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  /**
   * Destructor for the t8_search_with_batched_queries class.
   */
  virtual ~t8_search_with_batched_queries () = default;

 private:
  /**
   * Checks if the search should stop due to empty queries.
   */
  bool
  stop_due_to_queries () override
  {
    return this->active_queries.empty ();
  }

  /**
   * Checks the queries against the current element.
   * \param[in] new_active_queries A reference to a vector of new active query indices.
   * \param[in] ltreeid The ID of the local tree.
   * \param[in] element A pointer to the current element being checked.
   * \param[in] is_leaf A boolean indicating if the current element is a leaf.
   * \param[in] leaf_elements A pointer to an array of leaf elements.
   * \param[in] tree_leaf_index The index of the current element in the leaf array.
   */
  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements,
                 const t8_locidx_t tree_leaf_index) override
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!this->active_queries.empty ()) {
      std::vector<bool> query_matches (this->queries.size ());
      this->queries_callback (this->forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, this->queries,
                              this->active_queries, query_matches, this->user_data);
      std::copy_if (this->active_queries.begin (), this->active_queries.end (), std::back_inserter (new_active_queries),
                    [&] (size_t &query_index) { return query_matches[query_index]; });
    }
  }

  /**
   * Updates the queries for the search.
   * \param[in] old_query_indices A reference to a vector of queries to be processed.
   */
  void
  update_queries (std::vector<size_t> &old_query_indices) override
  {
    this->active_queries = old_query_indices;
  }

  /**
   * Initializes the queries for the search.
   */
  void
  init_queries () override
  {
    this->active_queries.resize (this->queries.size ());
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  }

  /**
   * The callback function for processing batched queries.
   */
  t8_search_batched_queries_callback<Query_T, Udata> queries_callback;
  /**
   * A reference to a vector of queries to be processed.
   */
  std::vector<Query_T> &queries;
  /**
   * A vector of indices of active queries.
   */
  std::vector<size_t> active_queries;
};

/**
 * A class that performs a search in the partition of a forest.
 */
class t8_partition_search_base {
 public:
  /**  Constructor for the t8_partition_search_base class.
   *
   *
   * This constructor initializes a t8_partition_search_base object with the
   * given forest. If the forest is not null, it increments the reference count
   * of the forest and asserts that the forest is committed.
   *
   * \param[in] forest A pointer to a t8_forest_t object. Defaults to nullptr.
   */
  t8_partition_search_base (t8_forest_t forest = nullptr): forest (forest)
  {
    if (forest != nullptr) {
      t8_forest_ref (forest);
      T8_ASSERT (t8_forest_is_committed (forest));
    }
  }

  /**  Update the forest for the search.
   *
   * This function updates the forest for the search. If the current forest is not null,
   * it decrements the reference count of the forest. It then asserts that the new forest
   * is not null and is committed. Finally, it increments the reference count of the new forest.
   *
   * \param[in] forest A pointer to a t8_forest_t object.
   */
  void
  update_forest (t8_forest_t forest)
  {
    if (this->forest != nullptr) {
      t8_forest_unref (&(this->forest));
    }
    T8_ASSERT (forest != nullptr);
    T8_ASSERT (t8_forest_is_committed (forest));
    t8_forest_ref (forest);
    this->forest = forest;
  }

  /**   Destructor for the t8_partition_search_base class.
   *
   * This destructor decrements the reference count of the forest if it is not null.
   */
  ~t8_partition_search_base ()
  {
    if (this->forest != nullptr) {
      t8_forest_unref (&(this->forest));
    }
  }

  /**   Perform the search.
   *
   * This function performs the search in the forest.
   */
  void
  do_search ();

  /**
   * A pointer to the forest whose partition is searched.
   */
  t8_forest_t forest;

 private:
  /**  Checks if the search should stop due to empty queries.
   *
   */
  virtual bool
  stop_due_to_queries ()
    = 0;

  /**  Checks an element during the search.
   *
   * This function is called for each element encountered during the search.
   * It passes the arguments to the callback function provided by the user.
   *
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   *
   * \return True if the search should continue, false otherwise.
   */
  virtual bool
  check_element (const t8_locidx_t ltreeid, const t8_element_t *element, const int pfirst, const int plast)
    = 0;

  /**  Checks queries during the search.
   *
   * This function is called to check queries during the search.
   * It passes the arguments to the callback function provided by the user.
   *
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   */
  virtual void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const int pfirst, const int plast)
    = 0;

  /**
   * Function the gives the user the opportunity to update the queries after each step in the recursion.
   *
   * \param old_query_indices
   */
  virtual void
  update_queries (std::vector<size_t> &old_query_indices)
    = 0;

  /**
   * Function gives the user the opportunity to set the queries to the initial full set before searching each tree.
   */
  virtual void
  init_queries ()
    = 0;
};

/**
 * A class that performs a search in the partition of a forest and supports the usage of user data.
 * \tparam Udata 
 */
template <typename Udata = void>
class t8_partition_search: public t8_partition_search_base {
 public:
  /**
  * Constructor for the t8_partition_search class.
  * This constructor initializes a t8_partition_search object with the given element callback,
  * forest, and user data. If the forest is not null, it increments the reference count
  * of the forest and asserts that the forest is committed.
  * 
  * \param[in] element_callback A callback function of type t8_partition_search_element_callback<Udata>.
  * \param[in] forest A pointer to a t8_forest_t object. Defaults to nullptr.
  * \param[in] user_data A pointer to a Udata object. Defaults to nullptr.
  */
  t8_partition_search (t8_partition_search_element_callback<Udata> element_callback, t8_forest_t forest = nullptr,
                       Udata *user_data = nullptr)
    : t8_partition_search_base (forest), element_callback (element_callback)
  {
    if (user_data != nullptr) {
      this->user_data = user_data;
    }
  }

  /**  Updates the user data associated with the object.
   *
   * This function sets the user data pointer to the provided Udata object.
   *
   * \param[in] udata A pointer to the Udata object to be set as the user data.
   */
  void
  update_user_data (Udata *udata)
  {
    if (udata != nullptr) {
      this->user_data = udata;
    }
  }

  /**
   * A pointer to the user data associated with the object.
   */
  Udata *user_data;

 private:
  /**
  *  In the partition search without queries, the search never stops due to empty queries.
  */
  bool
  stop_due_to_queries () override
  {
    return false;
  }

  /**
   *  In the partition search without queries, the element check simply forwards the call to the user callback.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   * \return True if the search should continue, false otherwise.
   */
  bool
  check_element (const t8_locidx_t ltreeid, const t8_element_t *element, const int pfirst, const int plast) override
  {
    T8_ASSERT (t8_forest_is_committed (this->forest));
    return this->element_callback (this->forest, ltreeid, element, pfirst, plast, user_data);
  }

  /**
   *  In the partition search without queries, no queries are checked.
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   */
  void
  check_queries ([[maybe_unused]] std::vector<size_t> &new_active_queries, [[maybe_unused]] const t8_locidx_t ltreeid,
                 [[maybe_unused]] const t8_element_t *element, [[maybe_unused]] const int pfirst,
                 [[maybe_unused]] const int plast) override
  {
    return;
  }

  /**
   *  In the partition search without queries, no queries are updated.
   * \param[in] old_query_indices A vector of indices of old queries.
   */
  void
  update_queries ([[maybe_unused]] std::vector<size_t> &old_query_indices)
  {
    return;
  }

  /**
   *  In the partition search without queries, no queries are initialized.
   */
  void
  init_queries () override
  {
    return;
  }

  /**
   *  The element callback function.
   */
  t8_partition_search_element_callback<Udata> element_callback;
};

/**
 *  A class that performs a search in the partition of a forest with queries.
 * Uses a filter-view to filter out the active queries. It is recommended to use this version of the search
 * if the number of queries is small or if the queries do not need any further computations to be evaluated.
 *
 * \tparam Query_T The type of queries
 * \tparam Udata The type of the user data, defaults to void.
*/
template <typename Query_T, typename Udata = void>
class t8_partition_search_with_queries: public t8_partition_search<Udata> {
 public:
  /**
   * Constructor for the t8_partition_search_with_queries class.
   * This constructor initializes a t8_partition_search_with_queries object with the given
   * element callback, queries callback, forest, and user data. If the forest is not null, it increments
   * the reference count of the forest and asserts that the forest is committed.
   * 
   * \param[in] element_callback A callback function of type t8_partition_search_element_callback<Udata>.
   * \param[in] queries_callback A callback function of type t8_partition_search_query_callback<Query_T, Udata>.
   * \param[in] queries A reference to a vector of Query_T objects.
   * \param[in] forest A pointer to a t8_forest_t object. Defaults to nullptr.
   * \param[in] user_data A pointer to a Udata object. Defaults to nullptr.
   */
  t8_partition_search_with_queries (t8_partition_search_element_callback<Udata> element_callback,
                                    t8_partition_search_query_callback<Query_T, Udata> queries_callback,
                                    std::vector<Query_T> &queries, const t8_forest_t forest = nullptr,
                                    Udata *user_data = nullptr)
    : t8_partition_search<Udata> (element_callback, forest, user_data), queries_callback (queries_callback),
      queries (queries)
  {
  }
  /**
   * Updates the queries associated with the object.
   *
   * \param[in] queries A reference to a vector of Query_T objects.
   */
  void
  update_queries (std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  /**
   * Destructor for the t8_partition_search_with_queries class.
   */
  ~t8_partition_search_with_queries ()
  {
  }

 private:
  /**
   * Stops the search if there are no active queries.
   *
   * \return True if the search should stop, false otherwise.
   */
  bool
  stop_due_to_queries () override
  {
    return this->active_queries.empty ();
  }

  /**
   * Checks queries during the search by invoking the user-provided callback function.
   *
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   */
  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const int pfirst, const int plast) override
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!this->active_queries.empty ()) {
      std::copy_if (this->active_queries.begin (), this->active_queries.end (), std::back_inserter (new_active_queries),
                    [&] (size_t &query_index) {
                      return this->queries_callback (this->forest, ltreeid, element, pfirst, plast,
                                                     queries[query_index], this->user_data);
                    });
    }
  }

  /**
   * Updates the queries associated with the object.
   *
   * \param[in] old_query_indices A vector of indices of old queries.
   */
  void
  update_queries (std::vector<size_t> &old_query_indices)
  {
    this->active_queries = old_query_indices;
  }

  /**
   * Initializes the queries associated with the object.
   */
  void
  init_queries () override
  {
    this->active_queries.resize (this->queries.size ());
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  }
  /**
   * Callback function for processing queries.
   */
  t8_partition_search_query_callback<Query_T, Udata> queries_callback;
  /**
   * A reference to a vector of queries.
   */
  std::vector<Query_T> &queries;
  /**
   * A vector of indices of active queries.
   */
  std::vector<size_t> active_queries;
};

/**
 * A class that performs a search in the partition of a forest with batched queries.
 *
 * All active queries are passed to the callback function, which processes them in a batch. It is recommended to
 * use this version of the searcch if further computations have to be done to evaluate the queries. That way these
 * precomputations are not done for every call to the callback again and only have to be evaluated once per call.
 *
 * \tparam Query_T The type of queries
 * \tparam Udata The type of the user data, defaults to void.
 */
template <typename Query_T, typename Udata = void>
class t8_partition_search_with_batched_queries: public t8_partition_search<Udata> {
 public:
  /**
  * Constructor for the t8_partition_search_with_batched_queries class.
  * This constructor initializes a t8_partition_search_with_batched_queries object with the given
  * element callback, queries callback, forest, and user data. If the forest is not null, it increments
  * the reference count of the forest and asserts that the forest is committed.
  * 
  * \param[in] element_callback A callback function of type t8_partition_search_element_callback<Udata>.
  * \param[in] queries_callback A callback function of type t8_partition_search_batched_queries_callback<Query_T, Udata>.
  * \param[in] queries A reference to a vector of Query_T objects.
  * \param[in] forest A pointer to a t8_forest_t object. Defaults to nullptr.
  * \param[in] user_data A pointer to a Udata object. Defaults to nullptr.
  */
  t8_partition_search_with_batched_queries (
    t8_partition_search_element_callback<Udata> element_callback,
    t8_partition_search_batched_queries_callback<Query_T, Udata> queries_callback, std::vector<Query_T> &queries,
    const t8_forest_t forest = nullptr, Udata *user_data = nullptr)
    : t8_partition_search<Udata> (element_callback, forest, user_data), queries_callback (queries_callback),
      queries (queries)
  {
  }
  /**
   * Updates the queries associated with the object.
   */
  void
  update_queries (std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  /**
   * Destructor for the t8_partition_search_with_batched_queries class.
   */
  ~t8_partition_search_with_batched_queries ()
  {
  }

 private:
  /**
   * In the partition search with batched queries, the search stops if there are no active queries left.
   * \return True if the search should stop, false otherwise.
   */
  bool
  stop_due_to_queries () override
  {
    return this->active_queries.empty ();
  }

  /**
   * Checks queries during the search by invoking the user-provided callback function.
   * \param[in] new_active_queries A vector of indices of active queries.
   * \param[in] ltreeid The local tree ID of the current element.
   * \param[in] element A pointer to the current element being processed.
   * \param[in] pfirst The first processor that owns part of \a element.
   * \param[in] plast The last processor that owns part of \a element.
   */
  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const int pfirst, const int plast) override
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!this->active_queries.empty ()) {
      std::vector<bool> query_matches (this->queries.size ());
      this->queries_callback (this->forest, ltreeid, element, pfirst, plast, this->queries, this->active_queries,
                              query_matches, this->user_data);
      std::copy_if (this->active_queries.begin (), this->active_queries.end (), std::back_inserter (new_active_queries),
                    [&] (size_t &query_index) { return query_matches[query_index]; });
    }
    std::swap (new_active_queries, this->active_queries);
  }

  /**
   * Updates the queries associated with the object.
   * \param[in] old_query_indices A vector of indices of old queries.
   */
  void
  update_queries (std::vector<size_t> &old_query_indices)
  {
    this->active_queries = old_query_indices;
  }
  /**
   * Initializes the queries associated with the object.
   */
  void
  init_queries () override
  {
    this->active_queries.resize (this->queries.size ());
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  }

  /**
   * Callback function for processing batched queries.
   */
  t8_partition_search_batched_queries_callback<Query_T, Udata> queries_callback;

  /**
   * A reference to a vector of queries.
   */
  std::vector<Query_T> &queries;

  /**
   * A vector of indices of active queries.
   */
  std::vector<size_t> active_queries;
};

#endif  // T8_FOREST_SEARCH_HXX
