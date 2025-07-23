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

/** \file t8_cmesh_cad.hxx
 * This geometry implements OpenCASCADE geometries. It enables the option to link different 
 * 1 and 2 dimensional cad geometries to the edges and faces of refinement trees. 
 * The geometry of the refinement tree is extended into the volume accordingly.
 */

#ifndef T8_CMESH_CAD
#define T8_CMESH_CAD

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <TopoDS.hxx>
#include <unordered_map>
#include <unordered_set>
#include <array>

struct t8_geom_data
{
  int entity_dim;
  int entity_tag;
  std::array<double, 2> location_on_curve;
};

class t8_boundary_node_geom_data_map {
 public:
  t8_boundary_node_geom_data_map (TopoDS_Shape& shape_in, t8_cmesh_t cmesh_in, double tolerance);

  std::unordered_map<t8_gloidx_t, t8_geom_data>
  get_boundary_node_geom_data_map ();

 private:
  void
  compute_geom_data_map ();

  TopoDS_Shape& shape;
  t8_cmesh_t cmesh;
  double tolerance;
  t8_geom_data geom_data;
  std::unordered_set<t8_gloidx_t> boundary_node_list;

  std::unordered_map<t8_gloidx_t, t8_geom_data> boundary_node_geom_data_map;
};

#endif /* T8_CMESH_CAD */
