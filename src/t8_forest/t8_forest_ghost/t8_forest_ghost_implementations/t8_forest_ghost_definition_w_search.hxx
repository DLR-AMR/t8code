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
#include <t8_forest/t8_forest_types.h>
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
  /** Base constructor with no arguments. We need this since it
   * is called from derived class constructors. */
  t8_forest_ghost_definition_w_search ()
  {
  }

  /**
   * Constructor with a search_function. Sets the type to T8_GHOST_TYPE_USER_DEFINED.
   * If do_ghost is called on this object,
   * the ghost layer will be created with a tree-based search (t8_forest_search)
   * with \a search_function as callback function.
   * \param search_function   The function used for the call bac
   */
  explicit t8_forest_ghost_definition_w_search (t8_forest_search_fn search_function,
                                                t8_forest_ghost_search_data *search_data)
    : t8_forest_ghost_definition (T8_GHOST_USER_DEFINED), search_fn (search_function), search_data (search_data)
  {
    T8_ASSERT (search_function != nullptr);
  }

  virtual ~t8_forest_ghost_definition_w_search ()
  {
    if (search_data != nullptr)
      delete search_data;
  }

  /** Create one layer of ghost elements for a forest.
   * \param [in,out]    forest     The forest.
   * \return T8_SUBROUTINE_SUCCESS if successful, T8_SUBROUTINE_FAILURE if not.
   * \a forest must be committed before calling this function.
   */
  virtual bool
  do_ghost (t8_forest_t forest) override;

 protected:
  /**
   * Fills the remote ghosts using a tree-based search.
   * \param [in,out]    forest     The forest.
   */
  virtual void
  search_for_ghost_elements (t8_forest_t forest);

  /**
   * Constructor for the derivided classes to set the type and the search_function.
   * \param [in] ghost_type       The type (faces, edges, user defined, ...) of the ghost_definition
   * \param [in] search_function  Function of type t8_forest_search_fn, used as callback function in search_for_ghost_elements
   * \param [in] search_data      Persistent data which can be used during the search. Ghost takes ownership of the data.
   */
  t8_forest_ghost_definition_w_search (const t8_ghost_type_t ghost_type, const t8_forest_search_fn search_function,
                                       t8_forest_ghost_search_data *search_data)
    : t8_forest_ghost_definition (ghost_type), search_fn (search_function), search_data (search_data)
  {
    T8_ASSERT (ghost_type != T8_GHOST_NONE);
  }

  t8_forest_search_fn search_fn {};         /** Callback function for t8_forest_search in search_for_ghost_elements */
  t8_forest_ghost_search_data *search_data; /** Persistent data which can be accessed during the search */
};

#endif /* !T8_FOREST_GHOST_DEFINITION_W_SEARCH_HXX */
