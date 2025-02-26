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

#ifndef T8_GHOST_DEFINITION_FACE_H
#define T8_GHOST_DEFINITION_FACE_H

#include <t8_forest/t8_forest_ghost_definition.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_iterate.h>  // Definition of t8_forest_search_fn

struct t8_forest_ghost_w_search: public t8_forest_ghost_definition
{
 public:
  /**
    * Constructors:
    * there are three ways to construct a object of t8_forest_ghost_w_search
    * t8_forest_ghost_w_search (), 
    * t8_forest_ghost_w_search (t8_forest_search_fn search_function),
    * t8_forest_ghost_w_search (const t8_ghost_type_t ghost_type)
    */
  t8_forest_ghost_w_search ();

  /**
   * Constructr of t8_forest_ghost_w_search by search_function
   * If do_ghost is called on this object, 
   * the ghost layer will be created by an treesearch (t8_forest_search)
   * with search_function as callbackfunction.
   * \note the t8_ghost_type_t of the object will we userdefined
   */
  explicit t8_forest_ghost_w_search (t8_forest_search_fn search_function)
    : t8_forest_ghost_definition (T8_GHOST_USERDEFINED), search_fn (search_function)
  {
    T8_ASSERT (search_function != nullptr);
  }

  /**
   * Constructr of t8_forest_ghost_w_search by type
   * The search_function is chosen by the type
   * \note currently only the type face is supported
   */
  explicit t8_forest_ghost_w_search (const t8_ghost_type_t ghost_type);

  virtual ~t8_forest_ghost_w_search ()
  {
  }

  /** Create one layer of ghost elements for a forest.
   * \param [in,out]    forest     The forest.
   * \a forest must be committed before calling this function.
   */
  virtual bool
  do_ghost (t8_forest_t forest) override;

 protected:
  /**
     * Equal to t8_forest_ghost_fill_remote_v3
     * so no support for version 1 and 2 of face neighbors.
     * Only the search_fn parameter for t8_forest_search 
     * is not the same as in t8_forest_ghost_fill_remote_v3.
     * Use the member variable of the class.
    */
  virtual void
  search_for_ghost_elements (t8_forest_t forest);

  /**
   * Constructor for the derivided classes to set the type and the search_function.
   * \param [in] ghost_type       The type (faces, edges, userdefind, ...) of the ghost_definition
   * \param [in] search_function  Function of type t8_forest_search_fn, used as callback function in search_for_ghost_elements
   */
  t8_forest_ghost_w_search (const t8_ghost_type_t ghost_type, const t8_forest_search_fn search_function)
    : t8_forest_ghost_definition (ghost_type), search_fn (search_function)
  {
    T8_ASSERT (ghost_type != T8_GHOST_NONE);
  }
  /** Callback function for t8_forest_search in search_for_ghost_elements */
  t8_forest_search_fn search_fn {};
};

struct t8_forest_ghost_face: public t8_forest_ghost_w_search
{
 public:
  /**
   * Constructor for the ghost class face.
   * do_ghost will construct a ghost layer with face neighbors
   * \param [in] version    one of tree versions (1,2,3) can be used
   * \note version 3 is the same treesearch as in t8_forest_ghost_w_search
   */
  explicit t8_forest_ghost_face (const int version);

  inline int
  get_version () const
  {
    return version;
  }

 protected:
  /**
   * Equal to t8_forest_ghost_fill_remote_v3 for version = 3
   * and t8_forest_ghost_fill_remote for version 1 and 2
   */
  void
  search_for_ghost_elements (t8_forest_t forest) override;

 private:
  int version {};
};

#endif /* !T8_GHOST_DEFINITION_FACE_H */
