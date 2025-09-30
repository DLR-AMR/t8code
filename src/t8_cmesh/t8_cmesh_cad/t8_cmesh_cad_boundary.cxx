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
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp.hxx>
#include <BRep_Tool.hxx>
#include <Standard_Real.hxx>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

t8_boundary_node_geom_data_map::t8_boundary_node_geom_data_map (TopoDS_Shape& shape_in, t8_cmesh_t cmesh_in,
                                                                double tolerance)
  : shape (shape_in), cmesh (cmesh_in), tolerance (tolerance)
{
  T8_ASSERT (cmesh->boundary_node_list != nullptr);
  TopExp::MapShapes (shape, TopAbs_VERTEX, cad_shape_vertex_map);
  TopExp::MapShapes (shape, TopAbs_EDGE, cad_shape_edge_map);
  TopExp::MapShapes (shape, TopAbs_FACE, cad_shape_face_map);
  boundary_node_list = cmesh->boundary_node_list->get_boundary_node_list ();
  T8_ASSERT (boundary_node_list.size () != 0);
  t8_debugf ("Boundary Node List Size = %ld\n", boundary_node_list.size ());
  compute_geom_data_map ();
}

void
t8_boundary_node_geom_data_map::compute_geom_data_map ()
{

  /* Create list of bounding boxes for curves */
  std::vector<Bnd_Box> edge_bboxes (cad_shape_edge_map.Extent () + 1);
  for (auto edge_iter = cad_shape_edge_map.cbegin (); edge_iter != cad_shape_edge_map.cend (); ++edge_iter) {
    const TopoDS_Edge& edge = TopoDS::Edge (*edge_iter);

    if (!BRep_Tool::Degenerated (edge)) {
      Bnd_Box box;
      BRepBndLib::Add (static_cast<const TopoDS_Shape&> (edge), box);

      int index = cad_shape_edge_map.FindIndex (*edge_iter);
      edge_bboxes[index] = box;
    }
  }

  /* Create list of bounding boxes for surfaces */
  std::vector<Bnd_Box> face_bboxes (cad_shape_face_map.Extent () + 1);
  for (auto face_iter = cad_shape_face_map.cbegin (); face_iter != cad_shape_face_map.cend (); ++face_iter) {
    const TopoDS_Face& face = TopoDS::Face (*face_iter);

    Bnd_Box box;
    BRepBndLib::Add (static_cast<const TopoDS_Shape&> (face), box);

    int index = cad_shape_face_map.FindIndex (*face_iter);
    face_bboxes[index] = box;
  }

  /* Iterate through t8_cmesh_boundary_node_list */
  for (auto bnl_iter = boundary_node_list.begin (); bnl_iter != boundary_node_list.end (); ++bnl_iter) {
    const tree_vertex_list tree_list = cmesh->vertex_connectivity->vertex_to_trees (*bnl_iter);
    t8_locidx_t local_tree_id = tree_list.at (0).first;
    int local_vertex_id = tree_list.at (0).second;
    double* vertices = (double*) t8_cmesh_get_tree_vertices (cmesh, local_tree_id);

    /* Get mesh node coordinates */
    const gp_Pnt mesh_pt (vertices[3 * local_vertex_id],      /* x-coordinate */
                          vertices[3 * local_vertex_id + 1],  /* y-coordinate */
                          vertices[3 * local_vertex_id + 2]); /* z-coordinate */

    /* Iterate through vertices of geometry */
    auto vertex_iter = cad_shape_vertex_map.cbegin ();
    for (; vertex_iter != cad_shape_vertex_map.cend (); ++vertex_iter) {
      int index = cad_shape_vertex_map.FindIndex (*vertex_iter);
      const gp_Pnt pt = BRep_Tool::Pnt (TopoDS::Vertex (*vertex_iter));
      if (mesh_pt.Distance (pt) <= tolerance) { /* If mesh node within tolerance of vertex */
        const t8_geom_data gd { 0, index, { -1, -1 } };
        boundary_node_geom_data_map.insert ({ *bnl_iter, gd }); /* append {global ID, t8_geom_data} to map */
        break;
      }
    }
    if (vertex_iter != cad_shape_vertex_map.cend ()) {
      continue;
    }

    /* Iterate through curves of geometry */
    auto edge_iter = cad_shape_edge_map.cbegin ();
    for (; edge_iter != cad_shape_edge_map.cend (); ++edge_iter) {
      int index = cad_shape_edge_map.FindIndex (*edge_iter);
      const TopoDS_Edge& edge = TopoDS::Edge (*edge_iter);

      /* Check if mesh node within bounding box */
      if (!BRep_Tool::Degenerated (edge) && !edge_bboxes[index].IsOut (mesh_pt)) {
        Standard_Real first, last;
        Handle (Geom_Curve) curve = BRep_Tool::Curve (edge, first, last);
        GeomAPI_ProjectPointOnCurve proj (mesh_pt, curve, first, last);

        /* Check if projection was successful and mesh node within tolerance of curve*/
        if (proj.NbPoints () && proj.LowerDistance () <= tolerance) {
          t8_geom_data gd;
          gd.entity_dim = 1;
          gd.entity_tag = index;
          gd.location_on_curve = { proj.LowerDistanceParameter (), -1 };

          /* append {global ID, t8_geom_data} to map */
          boundary_node_geom_data_map.insert ({ *bnl_iter, gd });
          break;
        }
      }
    }

    if (edge_iter != cad_shape_edge_map.cend ()) {
      continue;
    }

    /* Iterate through surfaces of geometry */
    auto face_iter = cad_shape_face_map.cbegin ();
    for (; face_iter != cad_shape_face_map.cend (); ++face_iter) {
      int index = cad_shape_face_map.FindIndex (*face_iter);
      if (!face_bboxes[index].IsOut (mesh_pt)) {
        const TopoDS_Face& face = TopoDS::Face (*face_iter);
        Handle (Geom_Surface) surface = BRep_Tool::Surface (face);
        GeomAPI_ProjectPointOnSurf proj (mesh_pt, surface);
        proj.Perform (mesh_pt);

        /* Check if projection was successful and mesh node within tolerance of curve*/
        if (proj.NbPoints () && proj.LowerDistance () <= tolerance) {
          double u, v;
          proj.LowerDistanceParameters (u, v);
          t8_geom_data gd;
          gd.entity_dim = 2;
          gd.entity_tag = index;
          gd.location_on_curve = { u, v };

          /* append {global ID, t8_geom_data} to map */
          boundary_node_geom_data_map.insert ({ *bnl_iter, gd });
          break;
        }
      }
    }
  }
}

std::unordered_map<t8_gloidx_t, t8_geom_data>
t8_boundary_node_geom_data_map::get_boundary_node_geom_data_map ()
{
  return this->boundary_node_geom_data_map;
};
