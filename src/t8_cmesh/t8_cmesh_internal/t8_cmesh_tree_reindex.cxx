/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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
#include <t8_cmesh_tree_reindex.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <array>

std::array<double, 24>
bbox_bounds_to_hex_vertices (const double bounds[6])
{
  const double x_min = bounds[0];
  const double x_max = bounds[1];
  const double y_min = bounds[2];
  const double y_max = bounds[3];
  const double z_min = bounds[4];
  const double z_max = bounds[5];

  return {
    x_min, y_min, z_min,  // vertex 0
    x_max, y_min, z_min,  // vertex 1
    x_min, y_max, z_min,  // vertex 2
    x_max, y_max, z_min,  // vertex 3

    x_min, y_min, z_max,  // vertex 4
    x_max, y_min, z_max,  // vertex 5
    x_min, y_max, z_max,  // vertex 6
    x_max, y_max, z_max   // vertex 7
  };
}

std::map<t8_gloidx_t, t8_gloidx_t>
t8_cmesh_reindex_tree (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  double bounding_box[6];
  t8_cmesh_get_local_bounding_box (cmesh, bounding_box);

  t8_cmesh_t bbox_cmesh;

  t8_cmesh_init (&bbox_cmesh);
  t8_cmesh_set_tree_class (bbox_cmesh, 0, T8_ECLASS_HEX);
  std::array<double, 24> vertices = bbox_bounds_to_hex_vertices (bounding_box);

  t8_cmesh_set_tree_vertices (bbox_cmesh, 0, vertices.data (), 8);

  t8_cmesh_commit (bbox_cmesh, comm);
};
