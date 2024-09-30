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

#ifndef T8_GHOST_INTERFACE_FACE_H
#define T8_GHOST_INTERFACE_FACE_H

#include <t8_forest/t8_forest_ghost_interface.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_iterate.h>  // Definition of t8_forest_search_query_fn

struct t8_forest_ghost_w_search: public t8_forest_ghost_interface
{
 public:
  /**
     * Constructor
    */
  t8_forest_ghost_w_search ();

  explicit t8_forest_ghost_w_search (t8_forest_search_query_fn search_function)
    : t8_forest_ghost_interface (T8_GHOST_USERDEFINED), search_fn (search_function)
  {
    T8_ASSERT (search_function != nullptr);
  }

  explicit t8_forest_ghost_w_search (t8_ghost_type_t ghost_type);

  virtual ~t8_forest_ghost_w_search ()
  {
  }

  virtual void
  do_ghost (t8_forest_t forest) override;

 protected:
  /**
     * Equal to t8_forest_ghost_fill_remote_v3
     * so no support for version 1 and 2 of face heigbors
     * Only the search_fn parameter for search it not the same as in t8_forest_ghost_fill_remote_v3
     * use the member variable of the class
    */
  virtual void
  search_for_ghost_elements (t8_forest_t forest);

  t8_forest_ghost_w_search (t8_ghost_type_t ghost_type, t8_forest_search_query_fn search_function)
    : t8_forest_ghost_interface (ghost_type), search_fn (search_function)
  {
    T8_ASSERT (ghost_type != T8_GHOST_NONE);
  }

  t8_forest_search_query_fn search_fn {};
};

struct t8_forest_ghost_face: public t8_forest_ghost_w_search
{
 public:
  explicit t8_forest_ghost_face (int version);

  inline int
  get_version () const
  {
    return version;
  }

 protected:
  void
  search_for_ghost_elements (t8_forest_t forest) override;

 private:
  int version {};
};

#endif /* !T8_GHOST_INTERFACE_FACE_H */
