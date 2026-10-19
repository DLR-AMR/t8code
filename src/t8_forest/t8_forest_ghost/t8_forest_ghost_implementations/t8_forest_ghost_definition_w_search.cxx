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
#include <t8_forest/t8_forest_types.h>

void
t8_forest_ghost_definition_w_search::fill_remote_ghosts (t8_forest_t forest)
{
  /* Store any internal data that may reside on the forest */
  void *const store_t8code_data = forest->t8code_data;
  /* Set the internal data for the search routine */
  forest->t8code_data = search_data;
  /* Loop over the trees of the forest */
  t8_forest_search (forest, search_fn, nullptr, nullptr);

  /* Reset the internal data from before the search */
  forest->t8code_data = store_t8code_data;
}

t8_forest_ghost_search_data *
t8_forest_ghost_get_search_data (const t8_forest_t forest)
{
  return (t8_forest_ghost_search_data *) forest->t8code_data;
}
