/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_ghost_definition_w_search.hxx
 * Implements a class to define ghosts via a tree search function.
 */

#ifndef T8_FOREST_GHOST_DEFINITION_W_SEARCH_HXX
#define T8_FOREST_GHOST_DEFINITION_W_SEARCH_HXX

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_base.hxx>
#include <t8_forest/t8_forest_iterate.h>

/** Base for the user search data used in t8_forest_ghost_definition_w_search  */
struct t8_forest_ghost_search_data
{
  /** Base constructor */
  t8_forest_ghost_search_data () {};

  /** Destructor */
  virtual ~t8_forest_ghost_search_data () {};
};

/**
 * Base class for all ghost definitions which use a tree-based search algorithm.
 */
struct t8_forest_ghost_definition_w_search: public t8_forest_ghost_definition
{
 public:
  /**
   * Constructor with a search_function.
   * If do_ghost is called on this object,
   * the ghost layer will be created with a tree-based search (t8_forest_search)
   * with \a search_function as callback function.
   * \param [in] search_function   The function used for the callback.
   * \param [in] search_data       Persistent data which can be used during the search. Ghost takes ownership of the data.
   * \note \a search_data is reachable from within \a search_function via
   * \ref t8_forest_ghost_get_search_data (forest).
   */
  explicit t8_forest_ghost_definition_w_search (t8_forest_search_fn search_function,
                                                t8_forest_ghost_search_data *search_data)
    : search_fn (search_function), search_data (search_data)
  {
    T8_ASSERT (search_function != nullptr);
  }

  ~t8_forest_ghost_definition_w_search () override
  {
    delete search_data;
  }

 protected:
  /** Base constructor with no arguments. We need this since it
   * can be called from derived class constructors. */
  t8_forest_ghost_definition_w_search ()
  {
  }

  /**
   * Fills the remote ghosts using a tree-based search.
   * \param [in,out]    forest     The forest.
   */
  void
  fill_remote_ghosts (t8_forest_t forest) override;

  t8_forest_search_fn search_fn {}; /**< Callback function for t8_forest_search in fill_remote_ghosts */
  /** Persistent data which can be accessed during the search.
   * \note Reachable from within \a search_fn via \ref t8_forest_ghost_get_search_data */
  t8_forest_ghost_search_data *search_data {};
};

/** Retrieve the search data of the ghost definition that is currently driving a search on \a forest.
 * Call this from within your \ref t8_forest_search_fn callback when your ghost definition derives from
 * \ref t8_forest_ghost_definition_w_search.
 * \param [in] forest  The forest passed to the search callback.
 * \return             The \a search_data of the ghost definition driving the current search.
 */
t8_forest_ghost_search_data *
t8_forest_ghost_get_search_data (const t8_forest_t forest);

#endif /* !T8_FOREST_GHOST_DEFINITION_W_SEARCH_HXX */
