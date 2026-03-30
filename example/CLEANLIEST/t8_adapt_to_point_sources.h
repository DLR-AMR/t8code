/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2026 the developers

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

#include <t8_forest/t8_forest_general.h> /* forest definition and basic interface. */

/** Struct wrapping a point source coordinate and a
 *  flag indicating whether it is in this partition.
 **/
struct t8_point_source_t
{
  double coordinates[3] = {};  //< The coordinates of the point source.
  int is_inside_partition
    = {};  //< Will be set to true if the point sources lies inside this process' parallel partition.
};

/**
 * Additional user data that we process during search.
 * For each element we count the number of point sources that it contains
 * and we count the total number of elements that we constructed during search.
 */
struct t8_adapt_to_point_sources_user_data_t
{
  std::vector<int> *points_per_element;  //< For each element the number of point sources inside it.
  t8_locidx_t num_elements_searched;     //< The total number of elements created.
};

t8_forest_t
t8_adapt_to_point_sources (t8_forest_t forest, sc_array_t *points);
