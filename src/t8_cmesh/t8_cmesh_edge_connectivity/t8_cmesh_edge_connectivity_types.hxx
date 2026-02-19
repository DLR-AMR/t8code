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

/* This file is part of the public interface, while
 * src/t8_cmesh/t8_cmesh_edge_connectivity/t8_cmesh_edge_connectivity.hxx 
 * is not. */

/** \file t8_cmesh_edge_connectivity_types.hxx
 * In this file we collect the data types required for cmesh 
 * global edge id storage.
 */

#ifndef T8_CMESH_EDGE_CONNECTIVITY_TYPES_HXX
#define T8_CMESH_EDGE_CONNECTIVITY_TYPES_HXX

/** Variable type for (tree_id, tree_edge_id) pair */
using tree_edge_pair = std::pair<t8_locidx_t, int>;

/** list of tree edge pairs, each global edge id maps to
  * such a list. */
using tree_edge_list = std::vector<tree_edge_pair>;

/** The internal storage data type used for storing the edge to tree data. */
using ett_storage_type = std::unordered_map<t8_gloidx_t, tree_edge_list>;

#endif /* !T8_CMESH_EDGE_CONNECTIVITY_TYPES_HXX */
