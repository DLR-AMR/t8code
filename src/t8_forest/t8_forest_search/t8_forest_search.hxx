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
#include <t8_forest/t8_forest.h>  // Ensure t8_forest_t is defined
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
using t8_search_element_callback = std::function<bool (
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, Udata *user_data)>;

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
  const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
  const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index, std::vector<Query_T> &queries,
  std::vector<size_t> &active_query_indices, std::vector<bool> &query_matches, Udata *user_data)>;

class t8_search_base {
 public:
  /**  \brief Constructor for the t8_search_base class.
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

  /**  \brief Update the forest for the search.
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

  /**  \brief Destructor for the t8_search_base class.
   *
   * This destructor decrements the reference count of the forest if it is not null.
   */
  ~t8_search_base ()
  {
    if (this->forest != nullptr) {
      t8_forest_unref (&(this->forest));
    }
  }

  /**  \brief Perform the search.
   *
   * This function performs the search in the forest. 
   */
  void
  do_search ();

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

  /** \brief Recursively searches the tree.
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
  search_recursion (const t8_locidx_t ltreeid, t8_element_t *element, const t8_eclass_scheme_c *ts,
                    t8_element_array_t *leaf_elements, const t8_locidx_t tree_lindex_of_first_leaf);

  /** \brief Checks if the search should stop due to empty queries.
   * 
   */
  virtual bool
  stop_due_to_queries ()
    = 0;

  /** \brief Checks an element during the search.
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

  /** \brief Checks queries during the search.
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
};

template <typename Udata = void>
class t8_search: public t8_search_base {
 public:
  t8_search (t8_search_element_callback<Udata> element_callback, t8_forest_t forest = nullptr,
             Udata *user_data = nullptr)
    : t8_search_base (forest), element_callback (element_callback), user_data (user_data)
  {
  }

  /** \brief Updates the user data associated with the object.
   *
   * This function sets the user data pointer to the provided Udata object.
   *
   * \param[in] udata A pointer to the Udata object to be set as the user data.
   */
  void
  update_user_data (Udata *udata)
  {
    this->user_data = udata;
  }

  Udata *user_data;

 private:
  bool
  stop_due_to_queries () override
  {
    return false;
  }

  bool
  check_element (const t8_locidx_t ltreeid, const t8_element_t *element, const bool is_leaf,
                 const t8_element_array_t *leaf_elements, const t8_locidx_t tree_leaf_index) override
  {
    T8_ASSERT (t8_forest_is_committed (this->forest));
    return this->element_callback (this->forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, user_data);
  }

  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements,
                 const t8_locidx_t tree_leaf_index) override
  {
    return;
  }

  t8_search_element_callback<Udata> element_callback;
};

template <typename Query_T, typename Udata = void>
class t8_search_with_queries: public t8_search<Udata> {
 public:
  t8_search_with_queries (t8_search_element_callback<Udata> element_callback,
                          t8_search_queries_callback<Query_T, Udata> queries_callback, std::vector<Query_T> &queries,
                          const t8_forest_t forest = nullptr, Udata *user_data = nullptr)
    : t8_search<Udata> (element_callback, forest, user_data), queries_callback (queries_callback), queries (queries)
  {
    this->active_queries.resize (queries.size ());
    std::iota (this->active_queries.begin (), this->active_queries.end (), 0);
  }

  void
  update_queries (std::vector<Query_T> &queries)
  {
    this->queries = queries;
  }

  ~t8_search_with_queries ()
  {
  }

 private:
  bool
  stop_due_to_queries () override
  {
    return this->active_queries.empty ();
  }

  void
  check_queries (std::vector<size_t> &new_active_queries, const t8_locidx_t ltreeid, const t8_element_t *element,
                 const bool is_leaf, const t8_element_array_t *leaf_elements,
                 const t8_locidx_t tree_leaf_index) override
  {
    T8_ASSERT (new_active_queries.empty ());
    if (!this->active_queries.empty ()) {
      std::vector<bool> query_matches (active_queries.size ());
      this->queries_callback (this->forest, ltreeid, element, is_leaf, leaf_elements, tree_leaf_index, this->queries,
                              this->active_queries, query_matches, this->user_data);
      if (!is_leaf) {
        for (size_t iactive : this->active_queries) {
          if (query_matches[iactive]) {
            new_active_queries.push_back (iactive);
          }
        }
      }
    }
  }

  t8_search_queries_callback<Query_T, Udata> queries_callback;
  std::vector<Query_T> &queries;
  std::vector<size_t> active_queries;
};

#endif  // T8_FOREST_SEARCH_HXX
