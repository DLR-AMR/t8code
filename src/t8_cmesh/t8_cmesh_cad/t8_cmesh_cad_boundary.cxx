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

/** \file t8_cmesh_cad.cxx
 * This geometry implements OpenCASCADE geometries. It enables the option to link different 
 * 1 and 2 dimensional cad geometries to the edges and faces of refinement trees. 
 * The geometry of the refinement tree is extended into the volume accordingly.
 */

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_cad_boundary.hxx>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_tree_to_vertex.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_vertex_to_tree.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <Standard_Real.hxx>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

t8_boundary_node_geom_data_map::t8_boundary_node_geom_data_map (TopoDS_Shape &shape_in, t8_cmesh_t cmesh_in,
                                                                double tolerance)
  : shape (shape_in), cmesh (cmesh_in), tolerance (tolerance)
{
  T8_ASSERT (cmesh->boundary_node_list != nullptr);
  boundary_node_list = cmesh->boundary_node_list->get_boundary_node_list ();
  T8_ASSERT (boundary_node_list.size () != 0);
  printf ("size = %ld\n", boundary_node_list.size ());
  compute_geom_data_map ();
}

void
t8_boundary_node_geom_data_map::compute_geom_data_map ()
{
  TopExp_Explorer dora;
  int tag_count = 1;
  for (dora.Init (shape, TopAbs_VERTEX); dora.More (); dora.Next ()) {
    const gp_Pnt point = BRep_Tool::Pnt (TopoDS::Vertex (dora.Current ()));
    const double cad_x_val = point.X ();
    const double cad_y_val = point.Y ();
    const double cad_z_val = point.Z ();

    for (auto iter = boundary_node_list.begin (); iter != boundary_node_list.end (); ++iter) {
      const tree_vertex_list tree_list = cmesh->vertex_connectivity->vertex_to_trees (*iter);
      t8_locidx_t local_tree_id = tree_list.at (0).first;
      int local_vertex_id = tree_list.at (0).second;
      double *vertices = (double *) t8_cmesh_get_tree_vertices (cmesh, local_tree_id);
      const double cmesh_x_coords = vertices[3 * local_vertex_id];
      const double cmesh_y_coords = vertices[3 * local_vertex_id + 1];
      const double cmesh_z_coords = vertices[3 * local_vertex_id + 2];

      const double dist = abs (sqrt (pow (cmesh_x_coords, 2) + pow (cmesh_y_coords, 2) + pow (cmesh_z_coords, 2))
                               - sqrt (pow (cad_x_val, 2) + pow (cad_y_val, 2) + pow (cad_z_val, 2)));

      if (dist <= tolerance) {
        t8_geom_data temp_geom_data;
        temp_geom_data.entity_dim = 0;
        temp_geom_data.entity_tag = tag_count;
        temp_geom_data.location_on_curve = { -1, -1 };

        boundary_node_geom_data_map.try_emplace (*iter, temp_geom_data);
        tag_count += 1;
      }
    }
  }
  t8_debugf ("Vertices added\n");

  tag_count = 1;
  for (dora.Init (shape, TopAbs_EDGE); dora.More (); dora.Next ()) {
    Standard_Real first, last;
    if (!BRep_Tool::Degenerated (TopoDS::Edge (dora.Current ()))) {
      Handle (Geom_Curve) geomCurve = BRep_Tool::Curve (TopoDS::Edge (dora.Current ()), first, last);
      for (auto iter = boundary_node_list.begin (); iter != boundary_node_list.end (); ++iter) {
        const tree_vertex_list tree_list = cmesh->vertex_connectivity->vertex_to_trees (*iter);
        t8_locidx_t local_tree_id = tree_list.at (0).first;
        int local_vertex_id = tree_list.at (0).second;
        double *vertices = (double *) t8_cmesh_get_tree_vertices (cmesh, local_tree_id);
        const double cmesh_x_coords = vertices[3 * local_vertex_id];
        const double cmesh_y_coords = vertices[3 * local_vertex_id + 1];
        const double cmesh_z_coords = vertices[3 * local_vertex_id + 2];

        const gp_Pnt vertex (cmesh_x_coords, cmesh_y_coords, cmesh_z_coords);

        GeomAPI_ProjectPointOnCurve projection (vertex, geomCurve);
        projection.Perform (vertex);
        if (projection.NbPoints ()) {
          double dist = projection.LowerDistance ();
          if (dist <= tolerance) {
            t8_geom_data temp_geom_data;
            temp_geom_data.entity_dim = 1;
            temp_geom_data.entity_tag = tag_count;
            temp_geom_data.location_on_curve = { projection.LowerDistanceParameter (), -1 };

            boundary_node_geom_data_map.try_emplace (*iter, temp_geom_data);
            tag_count += 1;
          }
        }
      }
    }
  }

  tag_count = 1;
  for (dora.Init (shape, TopAbs_FACE); dora.More (); dora.Next ()) {
    Handle (Geom_Surface) surfer = BRep_Tool::Surface (TopoDS::Face (dora.Current ()));
    for (auto iter = boundary_node_list.begin (); iter != boundary_node_list.end (); ++iter) {
      const tree_vertex_list tree_list = cmesh->vertex_connectivity->vertex_to_trees (*iter);
      t8_locidx_t local_tree_id = tree_list.at (0).first;
      int local_vertex_id = tree_list.at (0).second;
      double *vertices = (double *) t8_cmesh_get_tree_vertices (cmesh, local_tree_id);
      const double cmesh_x_coords = vertices[3 * local_vertex_id];
      const double cmesh_y_coords = vertices[3 * local_vertex_id + 1];
      const double cmesh_z_coords = vertices[3 * local_vertex_id + 2];

      const gp_Pnt vertex (cmesh_x_coords, cmesh_y_coords, cmesh_z_coords);

      GeomAPI_ProjectPointOnSurf projection (vertex, surfer);
      projection.Perform (vertex);
      if (projection.NbPoints ()) {
        double dist = projection.LowerDistance ();
        if (dist <= tolerance) {
          double u;
          double v;
          projection.LowerDistanceParameters (u, v);
          t8_geom_data temp_geom_data;
          temp_geom_data.entity_dim = 2;
          temp_geom_data.entity_tag = tag_count;
          temp_geom_data.location_on_curve = { u, v };

          boundary_node_geom_data_map.try_emplace (*iter, temp_geom_data);
          tag_count += 1;
        }
      }
    }
  }
  t8_debugf ("Faces added\n");
}

std::unordered_map<t8_gloidx_t, t8_geom_data>
t8_boundary_node_geom_data_map::get_boundary_node_geom_data_map ()
{
  return this->boundary_node_geom_data_map;
};
