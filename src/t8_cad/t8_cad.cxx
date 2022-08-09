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

#include <t8_cad/t8_cad.hxx>
#include <t8_cad/t8_cad.h>
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
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

/* *INDENT-OFF* */
t8_cad::t8_cad (const char *fileprefix)
{
  BRep_Builder        builder;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  BRepTools::Read (occ_shape, is, builder);
  is.close ();
  if (occ_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape \n");
  }
  TopExp::MapShapes (occ_shape, TopAbs_VERTEX, occ_shape_vertex_map);
  TopExp::MapShapes (occ_shape, TopAbs_EDGE, occ_shape_edge_map);
  TopExp::MapShapes (occ_shape, TopAbs_FACE, occ_shape_face_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_VERTEX, TopAbs_EDGE,
                                       occ_shape_vertex2edge_map);
  TopExp::MapShapesAndUniqueAncestors (occ_shape, TopAbs_EDGE, TopAbs_FACE,
                                       occ_shape_edge2face_map);
  BRepBndLib::AddOBB (occ_shape, occ_shape_bounding_box);
  gp_XYZ xyz = occ_shape_bounding_box.Center();
  t8_productionf("center: %f %f %f \n"
                 "size: %f %f %f", xyz.X(), xyz.Y(), xyz.Z(),
                 occ_shape_bounding_box.XHSize(), occ_shape_bounding_box.YHSize(), occ_shape_bounding_box.ZHSize());
}

t8_cad::t8_cad (const TopoDS_Shape occ_shape)
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
  BRepBndLib::AddOBB (occ_shape, occ_shape_bounding_box);
}

gp_Pnt
t8_cad::t8_cad_get_occ_point (const int index) const
{
  T8_ASSERT (index <= occ_shape_vertex_map.Size());
  return BRep_Tool::Pnt(TopoDS::Vertex(occ_shape_vertex_map.FindKey(index)));
}

Handle_Geom_Curve
t8_cad::t8_cad_get_occ_curve (const int index) const
{
  T8_ASSERT (index <= occ_shape_edge_map.Size());
  Standard_Real first, last;
  return BRep_Tool::Curve(TopoDS::Edge(occ_shape_edge_map.FindKey(index)), 
                          first, last);
}

Handle_Geom_Surface
t8_cad::t8_cad_get_occ_surface (const int index) const
{
  T8_ASSERT (index <= occ_shape_face_map.Size());
  return BRep_Tool::Surface(TopoDS::Face(occ_shape_face_map.FindKey(index)));
}

TopTools_IndexedMapOfShape 
t8_cad::t8_cad_get_occ_shape_vertex_map() const
{
  return occ_shape_vertex_map;
}

TopTools_IndexedMapOfShape 
t8_cad::t8_cad_get_occ_shape_edge_map() const
{
  return occ_shape_edge_map;
}

TopTools_IndexedMapOfShape
t8_cad::t8_cad_get_occ_shape_face_map() const
{
  return occ_shape_face_map;
}

int
t8_cad::t8_cad_get_common_edge (const int vertex1_index, 
                                          const int vertex2_index) const
{
  auto collection1 = occ_shape_vertex2edge_map.FindFromIndex(vertex1_index);
  auto collection2 = occ_shape_vertex2edge_map.FindFromIndex(vertex2_index);

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
t8_cad::t8_cad_get_common_face (const int edge1_index, 
                                          const int edge2_index) const
{
  auto collection1 = occ_shape_edge2face_map.FindFromIndex(edge1_index);
  auto collection2 = occ_shape_edge2face_map.FindFromIndex(edge2_index);

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
t8_cad::t8_cad_is_vertex_on_edge (const int vertex_index, 
                                            const int edge_index) const
{
  auto collection = occ_shape_vertex2edge_map.FindFromIndex(vertex_index);
  return collection.Contains(occ_shape_edge_map.FindKey(edge_index));
}

int
t8_cad::t8_cad_is_edge_on_face (const int edge_index, 
                                          const int face_index) const
{
  auto collection = occ_shape_edge2face_map.FindFromIndex(edge_index);
  return collection.Contains(occ_shape_face_map.FindKey(face_index));
}

int
t8_cad::t8_cad_is_vertex_on_face (const int vertex_index, 
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
t8_cad::t8_cad_get_parameter_of_vertex_on_edge(const int vertex_index, 
                                                         const int edge_index, 
                                                         double* edge_param) const
{
  T8_ASSERT(t8_cad::t8_cad_is_vertex_on_edge(vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Edge edge = TopoDS::Edge(occ_shape_edge_map.FindKey(edge_index));
  *edge_param = BRep_Tool::Parameter(vertex, edge);
}

void 
t8_cad::t8_cad_get_parameters_of_vertex_on_face(const int vertex_index, 
                                                          const int face_index, 
                                                          double* face_params) const
{
  T8_ASSERT(t8_cad::t8_cad_is_vertex_on_face(vertex_index, 
                                                       face_index));
  gp_Pnt2d uv;
  TopoDS_Vertex vertex = TopoDS::Vertex(occ_shape_vertex_map.FindKey(vertex_index));
  TopoDS_Face face = TopoDS::Face(occ_shape_face_map.FindKey(face_index));
  uv = BRep_Tool::Parameters(vertex, face);
  face_params[0] = uv.X();
  face_params[1] = uv.Y();
}

void 
t8_cad::t8_cad_edge_parameter_to_face_parameters(const int edge_index, 
                                                           const int face_index, 
                                                           const double edge_param, 
                                                           double* face_params) const
{
  T8_ASSERT(t8_cad::t8_cad_is_edge_on_face(edge_index, face_index));
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

int
t8_cad::t8_cad_is_element_inside_shape(t8_forest_t forest,
                                       t8_locidx_t ltreeid, 
                                       const t8_element_t *element) const
{
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_eclass_t         tree_class, element_class;
  t8_eclass_scheme_c *ts;
  tree_class = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, tree_class);
  /* Check if element is valid */
  T8_ASSERT (ts->t8_element_is_valid (element));
  
  /* Check if element is quad or hex */
  element_class = ts->t8_element_shape(element);
  T8_ASSERT (element_class == T8_ECLASS_HEX || element_class == T8_ECLASS_QUAD);
  
#if T8_ENABLE_DEBUG
  /* Check if element is axis oriented */
  double corner_values[24];
  int num_equal_coordinates;
  for (int corner = 0; corner < t8_eclass_num_vertices[element_class]; ++corner) {
    ts->t8_element_vertex_reference_coords(element, 
                                           corner, 
                                           corner_values + corner * 3);
  }
  /* An element is axis oriented if all edges align to at least one axis */
  for (int edge = 0; edge < t8_eclass_num_edges[element_class]; ++edge) {
    num_equal_coordinates = 0;
    for (int dim = 0; dim < 3; ++dim) {
      if (element_class == T8_ECLASS_HEX) {
        if (std::abs(corner_values[t8_edge_vertex_to_tree_vertex[edge][0] * 3 + dim]
                     - corner_values[t8_edge_vertex_to_tree_vertex[edge][1] * 3 + dim])
            <= DBL_EPSILON) {
          ++num_equal_coordinates;
        }
      }
      else {
        if (std::abs(corner_values[t8_face_vertex_to_tree_vertex[T8_ECLASS_QUAD][edge][0] * 3 + dim]
                     - corner_values[t8_face_vertex_to_tree_vertex[T8_ECLASS_QUAD][edge][1] * 3 + dim])
            <= DBL_EPSILON) {
          ++num_equal_coordinates;
        }
      }
    }
    T8_ASSERT(num_equal_coordinates >= 2);
  }
#endif /* T8_ENABLE_DEBUG */
  /* Compute bounding box of element */
  double corner_coords[3];
  const int max_corner_number = t8_eclass_num_vertices[element_class];
  ts->t8_element_vertex_reference_coords(element, 0, corner_coords);
  gp_Pnt box_min = gp_Pnt(corner_coords[0], corner_coords[1], corner_coords[2]);
  ts->t8_element_vertex_reference_coords(element, max_corner_number, 
                                         corner_coords);
  gp_Pnt box_max = gp_Pnt(corner_coords[0], corner_coords[1], corner_coords[2]);
  Bnd_Box unoriented_bounding_box = Bnd_Box(box_min, box_max);
  Bnd_OBB element_bounding_box = Bnd_OBB(unoriented_bounding_box);
  element_bounding_box.SetAABox(1);
  
  /* Check if element bounding box is outside of shape bounding box. 
   * If true, element is completely outside of the shape. */
  if (occ_shape_bounding_box.IsOut(element_bounding_box)) {
    return 0;
  }

  return 1;

}

int
t8_cad::t8_cad_is_point_inside_shape (const double *coords, double tol) const
{
  gp_Pnt pnt = gp_Pnt(coords[0], coords[1], coords[2]);
  if (occ_shape_bounding_box.IsOut(pnt)) {
    return 0;
  }
  BRepClass3d_SolidClassifier classifier = BRepClass3d_SolidClassifier(occ_shape);
  classifier.Perform(pnt, tol);
  return classifier.State() ? 0 : 1;
}
/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
