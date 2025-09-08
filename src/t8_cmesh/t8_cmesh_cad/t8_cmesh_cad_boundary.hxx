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

/** \file t8_cmesh_cad_boundary.hxx
 * Enables the correct identification and remapping of boundary nodes on a geometry
 */

#ifndef T8_CMESH_CAD
#define T8_CMESH_CAD

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <TopoDS.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <unordered_map>
#include <unordered_set>
#include <array>

/**
 *  Struct to store geometric data
 */
struct t8_geom_data
{
  int entity_dim; /**< The entity's dimension (vertex->0, curve->1, surface->2) */
  int entity_tag; /**< The entity's key in its corresponding dimensional entity map */
  std::array<double, 2>
    location_on_curve; /**< Parameters if point was projected -> (0D -> {-1, -1}, 1D -> {u, -1}, 2D -> {u, v}) */
};

class t8_boundary_node_geom_data_map {
 public:
  /** Constructor from a given CAD geometry and a corresponding committed cmesh. 
  * 
  * \param [in] shape_in    A CAD geometry
  * \param [in] cmesh_in    A committed cmesh
  * \param [in] tolerance   A user defined tolerance to specify the precision of the boundary node remapping. Defaulted to 1e-7.
  * 
  */
  t8_boundary_node_geom_data_map (TopoDS_Shape& shape_in, t8_cmesh_t cmesh_in, const double tolerance = 1e-7);

  /** Getter function for the geometry data map
   * 
   * \return    A std::unordered_map with multiple pairs of global node ID as the key and a t8_geom_data struct as the value
   * 
   */
  std::unordered_map<t8_gloidx_t, t8_geom_data>
  get_boundary_node_geom_data_map ();

 private:
  /** Function to fill the \a t8_boundary_node_geom_data_map */
  void
  compute_geom_data_map ();

  TopoDS_Shape& shape;                                /** CAD Geometry */
  t8_cmesh_t cmesh;                                   /** Corresponding mesh */
  const double tolerance;                             /** User specified tolerance */
  std::unordered_set<t8_gloidx_t> boundary_node_list; /** Boundary node list of the given cmesh */
  std::unordered_map<t8_gloidx_t, t8_geom_data>
    boundary_node_geom_data_map; /** Hashmap with key-value pairs of global_index - geom_data */

  TopTools_IndexedMapOfShape cad_shape_vertex_map; /**< Map of all TopoDS_Vertex in shape. */
  TopTools_IndexedMapOfShape cad_shape_edge_map;   /**< Map of all TopoDS_Edge in shape. */
  TopTools_IndexedMapOfShape cad_shape_face_map;   /**< Map of all TopoDS_Face in shape. */
};

#endif /* T8_CMESH_CAD */
