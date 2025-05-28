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

/** \file t8_forest_ghost_definition_w_search.cxx
 * Implementations for t8_forest_ghost_definition_w_search.hxx
 */

#include <t8_forest/t8_forest_ghost/t8_forest_ghost_implementations/t8_forest_ghost_definition_w_search.hxx>
#include <t8_forest/t8_forest_ghost/t8_forest_ghost_definition_helpers.hxx>

bool
t8_forest_ghost_definition_w_search::do_ghost (t8_forest_t forest)
{

  if (t8_ghost_get_type () == T8_GHOST_NONE) {
    t8_debugf ("WARNING: Trying to construct ghosts with ghost_type NONE. "
               "Ghost layer is not constructed.\n");
    return T8_SUBROUTINE_FAILURE;
  }

  communicate_ownerships (forest);

  if (t8_forest_get_local_num_elements (forest) > 0) {

    /* Initialize the ghost structure */
    t8_forest_ghost_init (&forest->ghosts, ghost_type);

    search_for_ghost_elements (forest);

    communicate_ghost_elements (forest);
  }
  clean_up (forest);

  return T8_SUBROUTINE_SUCCESS;
}

void
t8_forest_ghost_definition_w_search::search_for_ghost_elements (t8_forest_t forest)
{
  void *store_user_data = NULL;

  /* Store any user data that may reside on the forest */
  store_user_data = t8_forest_get_user_data (forest);
  /* Set the user data for the search routine */
  t8_forest_set_user_data (forest, &search_data);
  /* Loop over the trees of the forest */
  t8_forest_search (forest, search_fn, NULL, NULL);

  /* Reset the user data from before search */
  t8_forest_set_user_data (forest, store_user_data);
}
