/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <t8_cad/t8_cad_geom.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
/* TODO: Remove t8_geometry_occ.hxx as soon as edge connectivity is defined in t8_eclass.h */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_occ.hxx>

#if T8_WITH_OCC
#include <BRep_Builder.hxx>
#include <TopExp.hxx>
#include <BRepTools.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

/* *INDENT-OFF* */
t8_cad_geom::t8_cad_geom (const char *fileprefix)
{
  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream filestream (current_file + ".brep");
  BRepTools::Read (occ_shape, filestream, builder);
  filestream.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape \n");
  }
  t8_cad_geom::t8_cad_init_internal_data ();
}

t8_cad_geom::t8_cad_geom (const TopoDS_Shape shape)
{
  occ_shape = shape;
  t8_cad_geom::t8_cad_init_internal_data ();
}

void
t8_cad_geom::t8_cad_init (const char *fileprefix)
{
  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream filestream (current_file + ".brep");
  BRepTools::Read (occ_shape, filestream, builder);
  filestream.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape \n");
  }
  t8_cad_geom::t8_cad_init_internal_data ();
}

void
t8_cad_geom::t8_cad_init (const TopoDS_Shape shape)
{
  occ_shape = shape;
  t8_cad_geom::t8_cad_init_internal_data ();
}

void
t8_cad_geom::t8_cad_init_internal_data ()
{
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Shape is null. \n");
  }
  TopExp::MapShapes (occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes (occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes (occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_VERTEX, TopAbs_EDGE,
                                       occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_EDGE, TopAbs_FACE,
                                       occ_shape_edge2face_map);
}

const gp_Pnt
t8_cad_geom::t8_cad_get_occ_point (const int index) const
{
  T8_ASSERT (index <= occ_shape_vertex_map.Size());
  return BRep_Tool::Pnt(TopoDS::Vertex(occ_shape_vertex_map.FindKey(index)));
}

const Handle_Geom_Curve
t8_cad_geom::t8_cad_get_occ_curve (const int index) const
{
  T8_ASSERT (index <= occ_shape_edge_map.Size());
  Standard_Real first, last;
  return BRep_Tool::Curve(TopoDS::Edge(occ_shape_edge_map.FindKey(index)), 
                          first, last);
}

const Handle_Geom_Surface
t8_cad_geom::t8_cad_get_occ_surface (const int index) const
{
  T8_ASSERT (index <= occ_shape_face_map.Size());
  return BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(index)));
}

const TopTools_IndexedMapOfShape 
t8_cad_geom::t8_cad_get_occ_shape_vertex_map() const
{
  return occ_shape_vertex_map;
}

const TopTools_IndexedMapOfShape 
t8_cad_geom::t8_cad_get_occ_shape_edge_map() const
{
  return occ_shape_edge_map;
}

const TopTools_IndexedMapOfShape
t8_cad_geom::t8_cad_get_occ_shape_face_map() const
{
  return occ_shape_face_map;
}

int
t8_cad_geom::t8_cad_get_common_edge (const int vertex1_index, 
                                     const int vertex2_index) const
{
  /* Save all edges the vertices lie on in a seperate collection */
  auto collection1 = occ_shape_vertex2edge_map.FindFromIndex(vertex1_index);
  auto collection2 = occ_shape_vertex2edge_map.FindFromIndex(vertex2_index);

  /* Iterate over each edge to check, if both vertices share an edge */
  for (auto edge1 = collection1.begin(); edge1 != collection1.end(); ++edge1)
  {
    for (auto edge2 = collection2.begin(); edge2 != collection2.end(); ++edge2)
    {
      if (edge1->IsEqual(*edge2))
      {
        return occ_shape_edge2face_map.FindIndex(*edge1);
      }
    }
  }
  return 0;
}

int
t8_cad_geom::t8_cad_get_common_face (const int edge1_index, 
                                     const int edge2_index) const
{
  /* Save all faces the edges lie on in a seperate collection */
  auto collection1 = occ_shape_edge2face_map.FindFromIndex(edge1_index);
  auto collection2 = occ_shape_edge2face_map.FindFromIndex(edge2_index);

  /* Iterate over each face to check, if both edges share a face */
  for (auto face1 = collection1.begin(); face1 != collection1.end(); ++face1)
  {
    for (auto face2 = collection2.begin(); face2 != collection2.end(); ++face2)
    {
      if (face1->IsEqual(*face2))
      {
        return occ_shape_face_map.FindIndex(*face1);
      }
    }
  }
  return 0;
}

int
t8_cad_geom::t8_cad_is_vertex_on_edge (const int vertex_index, 
                                       const int edge_index) const
{
  auto collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  return collection.Contains(occ_shape_edge_map.FindKey(edge_index));
}

int
t8_cad_geom::t8_cad_is_edge_on_face (const int edge_index, 
                                     const int face_index) const
{
  auto collection = occ_shape_edge2face_map.FindFromIndex(edge_index);
  return collection.Contains(occ_shape_face_map.FindKey(face_index));
}

int
t8_cad_geom::t8_cad_is_vertex_on_face (const int vertex_index, 
                                       const int face_index) const
{
  auto edge_collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  for (auto edge = edge_collection.begin(); edge != edge_collection.end(); ++edge)
  {
    auto face_collection = occ_shape_edge2face_map.FindFromKey(*edge);
    if (face_collection.Contains(occ_shape_face_map.FindKey(face_index)))
    {
      return 1;
    }
  }
  return 0;
}

void 
t8_cad_geom::t8_cad_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                                    const int edge_index, 
                                                    double* edge_param) const
{
  T8_ASSERT(t8_cad_geom::t8_cad_is_vertex_on_edge(vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  *edge_param = BRep_Tool::Parameter(vertex, edge);
}

void 
t8_cad_geom::t8_cad_get_parameters_of_vertex_on_face(const int vertex_index, 
                                                     const int face_index, 
                                                     double* face_params) const
{
  T8_ASSERT(t8_cad_geom::t8_cad_is_vertex_on_face(vertex_index, 
                                                       face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  uv = BRep_Tool::Parameters(vertex, face);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}

void 
t8_cad_geom::t8_cad_edge_parameter_to_face_parameters(const int edge_index, 
                                                      const int face_index, 
                                                      const double edge_param, 
                                                      double* face_params) const
{
  T8_ASSERT(t8_cad_geom::t8_cad_is_edge_on_face(edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  Handle_Geom2d_Curve curve_on_surface  = BRep_Tool::CurveOnSurface(edge, face, 
                                                                    first, last);
  curve_on_surface->D0(edge_param, uv);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
