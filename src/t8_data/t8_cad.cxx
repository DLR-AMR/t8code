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

#include <t8.h>
#include <t8_data/t8_cad.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#include <TopoDS.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <Standard_Version.hxx>

t8_cad::t8_cad (std::string fileprefix)
{
  BRep_Builder builder;
  std::ifstream is (fileprefix + ".brep");
  if (is.is_open () == false) {
    SC_ABORTF ("Cannot find the file %s.brep.\n", fileprefix.c_str ());
  }
  BRepTools::Read (cad_shape, is, builder);
  is.close ();
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape. "
               "The cad file may be written with a newer cad version. "
               "Linked cad version: %s",
               OCC_VERSION_COMPLETE);
  }
  TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
  TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
  TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_VERTEX, TopAbs_EDGE, cad_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_EDGE, TopAbs_FACE, cad_shape_edge2face_map);
}

t8_cad::t8_cad (const TopoDS_Shape cad_shape)
{
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
  TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
  TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_VERTEX, TopAbs_EDGE, cad_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_EDGE, TopAbs_FACE, cad_shape_edge2face_map);
}

t8_cad::t8_cad ()
{
  cad_shape.Nullify ();
}

int
t8_cad::t8_geom_is_line (const int curve_index) const
{
  const Handle_Geom_Curve curve = t8_geom_get_cad_curve (curve_index);
  const GeomAdaptor_Curve curve_adaptor (curve);
  return curve_adaptor.GetType () == GeomAbs_Line;
}

int
t8_cad::t8_geom_is_plane (const int surface_index) const
{
  const Handle_Geom_Surface surface = t8_geom_get_cad_surface (surface_index);
  const GeomAdaptor_Surface surface_adaptor (surface);
  return surface_adaptor.GetType () == GeomAbs_Plane;
}

const gp_Pnt
t8_cad::t8_geom_get_cad_point (const int index) const
{
  T8_ASSERT (index <= cad_shape_vertex_map.Size ());
  return BRep_Tool::Pnt (TopoDS::Vertex (cad_shape_vertex_map.FindKey (index)));
}

const Handle_Geom_Curve
t8_cad::t8_geom_get_cad_curve (const int index) const
{
  T8_ASSERT (index <= cad_shape_edge_map.Size ());
  Standard_Real first, last;
  return BRep_Tool::Curve (TopoDS::Edge (cad_shape_edge_map.FindKey (index)), first, last);
}

const Handle_Geom_Surface
t8_cad::t8_geom_get_cad_surface (const int index) const
{
  T8_ASSERT (index <= cad_shape_face_map.Size ());
  return BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (index)));
}

const TopTools_IndexedMapOfShape
t8_cad::t8_geom_get_cad_shape_vertex_map () const
{
  return cad_shape_vertex_map;
}

const TopTools_IndexedMapOfShape
t8_cad::t8_geom_get_cad_shape_edge_map () const
{
  return cad_shape_edge_map;
}

const TopTools_IndexedMapOfShape
t8_cad::t8_geom_get_cad_shape_face_map () const
{
  return cad_shape_face_map;
}

int
t8_cad::t8_geom_get_common_edge (const int vertex1_index, const int vertex2_index) const
{
  const TopTools_ListOfShape collection1 = cad_shape_vertex2edge_map.FindFromIndex (vertex1_index);
  const TopTools_ListOfShape collection2 = cad_shape_vertex2edge_map.FindFromIndex (vertex2_index);

  for (auto edge1 = collection1.begin (); edge1 != collection1.end (); ++edge1) {
    for (auto edge2 = collection2.begin (); edge2 != collection2.end (); ++edge2) {
      if (edge1->IsEqual (*edge2)) {
        return cad_shape_edge2face_map.FindIndex (*edge1);
      }
    }
  }
  return 0;
}

int
t8_cad::t8_geom_get_common_face (const int edge1_index, const int edge2_index) const
{
  const TopTools_ListOfShape collection1 = cad_shape_edge2face_map.FindFromIndex (edge1_index);
  const TopTools_ListOfShape collection2 = cad_shape_edge2face_map.FindFromIndex (edge2_index);

  for (auto face1 = collection1.begin (); face1 != collection1.end (); ++face1) {
    for (auto face2 = collection2.begin (); face2 != collection2.end (); ++face2) {
      if (face1->IsEqual (*face2)) {
        return cad_shape_face_map.FindIndex (*face1);
      }
    }
  }
  return 0;
}

int
t8_cad::t8_geom_is_vertex_on_edge (const int vertex_index, const int edge_index) const
{
  const TopTools_ListOfShape collection = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
  return collection.Contains (cad_shape_edge_map.FindKey (edge_index));
}

int
t8_cad::t8_geom_is_edge_on_face (const int edge_index, const int face_index) const
{
  const TopTools_ListOfShape collection = cad_shape_edge2face_map.FindFromIndex (edge_index);
  return collection.Contains (cad_shape_face_map.FindKey (face_index));
}

int
t8_cad::t8_geom_is_vertex_on_face (const int vertex_index, const int face_index) const
{
  const TopTools_ListOfShape edge_collection = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
  for (auto edge = edge_collection.begin (); edge != edge_collection.end (); ++edge) {
    const TopTools_ListOfShape face_collection = cad_shape_edge2face_map.FindFromKey (*edge);
    if (face_collection.Contains (cad_shape_face_map.FindKey (face_index))) {
      return 1;
    }
  }
  return 0;
}

void
t8_cad::t8_geom_get_parameter_of_vertex_on_edge (const int vertex_index, const int edge_index, double *edge_param) const
{
  T8_ASSERT (t8_cad::t8_geom_is_vertex_on_edge (vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));
  *edge_param = BRep_Tool::Parameter (vertex, edge);
}

void
t8_cad::t8_geom_get_parameters_of_vertex_on_face (const int vertex_index, const int face_index,
                                                  double *face_params) const
{
  T8_ASSERT (t8_cad::t8_geom_is_vertex_on_face (vertex_index, face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));
  uv = BRep_Tool::Parameters (vertex, face);
  face_params[0] = uv.X ();
  face_params[1] = uv.Y ();
}

void
t8_cad::t8_geom_edge_parameter_to_face_parameters (const int edge_index, const int face_index, const int num_face_nodes,
                                                   const double edge_param, const double *surface_params,
                                                   double *face_params) const
{
  T8_ASSERT (t8_cad::t8_geom_is_edge_on_face (edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));
  Handle_Geom2d_Curve curve_on_surface = BRep_Tool::CurveOnSurface (edge, face, first, last);
  Handle_Geom_Surface surface = BRep_Tool::Surface (face);
  curve_on_surface->D0 (edge_param, uv);
  face_params[0] = uv.X ();
  face_params[1] = uv.Y ();

  /* Check for right conversion of edge to surface parameter and correct if needed */
  /* Checking u parameter */
  if (surface_params != NULL) {
    double parametric_bounds[4];
    surface->Bounds (parametric_bounds[0], parametric_bounds[1], parametric_bounds[2], parametric_bounds[3]);
    if (surface->IsUClosed ()) {
      for (int i_face_node = 0; i_face_node < num_face_nodes; ++i_face_node) {
        if (surface_params[i_face_node * 2] == parametric_bounds[0]) {
          if (face_params[0] == parametric_bounds[1]) {
            face_params[0] = parametric_bounds[0];
          }
        }
        else if (surface_params[i_face_node * 2] == parametric_bounds[1]) {
          if (face_params[0] == parametric_bounds[0]) {
            face_params[0] = parametric_bounds[1];
          }
        }
      }
    }
    /* Checking v parameter */
    if (surface->IsVClosed ()) {
      for (int i_face_node = 0; i_face_node < num_face_nodes; ++i_face_node) {
        if (surface_params[i_face_node * 2 + 1] == parametric_bounds[0]) {
          if (face_params[1] == parametric_bounds[1]) {
            face_params[1] = parametric_bounds[0];
          }
        }
        else if (surface_params[i_face_node * 2 + 1] == parametric_bounds[1]) {
          if (face_params[1] == parametric_bounds[0]) {
            face_params[1] = parametric_bounds[1];
          }
        }
      }
    }
  }
}

void
t8_cad::t8_geom_get_face_parametric_bounds (const int surface_index, double *bounds) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (surface_index);
  cad_surface->Bounds (bounds[0], bounds[1], bounds[2], bounds[3]);
}

void
t8_cad::t8_geom_get_edge_parametric_bounds (const int edge_index, double *bounds) const
{
  const Handle_Geom_Curve cad_edge = t8_geom_get_cad_curve (edge_index);
  bounds[0] = cad_edge->FirstParameter ();
  bounds[1] = cad_edge->LastParameter ();
}

int
t8_cad::t8_geom_is_edge_closed (int edge_index) const
{
  const Handle_Geom_Curve cad_edge = t8_geom_get_cad_curve (edge_index);
  return cad_edge->IsClosed ();
}

int
t8_cad::t8_geom_is_surface_closed (int geometry_index, int parameter) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (geometry_index);
  switch (parameter) {
  case 0:
    return cad_surface->IsUClosed ();
    break;
  case 1:
    return cad_surface->IsVClosed ();
    break;
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
}
