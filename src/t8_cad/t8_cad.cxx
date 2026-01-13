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
#include <t8_cad/t8_cad.hxx>
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
#include <ShapeAnalysis_Edge.hxx>
#include <TopExp_Explorer.hxx>

#if T8_ENABLE_DEBUG
#include <Precision.hxx>
#endif

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
  TopExp::MapShapesAndUniqueAncestors (cad_shape, TopAbs_VERTEX, TopAbs_FACE, cad_shape_vertex2face_map);
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

const TopoDS_Vertex
t8_cad::t8_geom_get_cad_vertex (const int index) const
{
  T8_ASSERT (index <= cad_shape_vertex_map.Size ());
  return TopoDS::Vertex (cad_shape_vertex_map.FindKey (index));
}

const TopoDS_Edge
t8_cad::t8_geom_get_cad_edge (const int index) const
{
  T8_ASSERT (index <= cad_shape_edge_map.Size ());
  return TopoDS::Edge (cad_shape_edge_map.FindKey (index));
}

const TopoDS_Face
t8_cad::t8_geom_get_cad_face (const int index) const
{
  T8_ASSERT (index <= cad_shape_face_map.Size ());
  return TopoDS::Face (cad_shape_face_map.FindKey (index));
}

const gp_Pnt
t8_cad::t8_geom_get_cad_point (const int index) const
{
  return BRep_Tool::Pnt (t8_cad::t8_geom_get_cad_vertex (index));
}

const Handle_Geom_Curve
t8_cad::t8_geom_get_cad_curve (const int index) const
{
  Standard_Real first, last;
  return BRep_Tool::Curve (t8_cad::t8_geom_get_cad_edge (index), first, last);
}

const Handle_Geom_Surface
t8_cad::t8_geom_get_cad_surface (const int index) const
{
  return BRep_Tool::Surface (t8_cad::t8_geom_get_cad_face (index));
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
t8_cad::t8_geom_get_common_edge_of_vertices (const int vertex1_index, const int vertex2_index) const
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
t8_cad::t8_geom_get_common_face_of_edges (const int edge1_index, const int edge2_index) const
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
t8_cad::t8_geom_get_common_face_of_vertex_and_edge (const int vertex_index, const int edge_index) const
{
  const TopTools_ListOfShape edge_collection = cad_shape_edge2face_map.FindFromIndex (edge_index);
  for (auto face = edge_collection.begin (); face != edge_collection.end (); ++face) {
    const size_t face_index = cad_shape_face_map.FindIndex (*face);
    if (t8_geom_is_vertex_on_face (vertex_index, face_index)) {
      return face_index;
    }
  }
  return 0;
}

int
t8_cad::t8_geom_get_common_face_of_vertices (const int vertex1_index, const int vertex2_index) const
{
  const TopTools_ListOfShape collection1 = cad_shape_vertex2face_map.FindFromIndex (vertex1_index);
  const TopTools_ListOfShape collection2 = cad_shape_vertex2face_map.FindFromIndex (vertex2_index);

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

bool
t8_cad::t8_geom_vertex_is_on_seam (const int vertex_index, const int face_index) const
{
  const auto face = t8_cad::t8_geom_get_cad_face (face_index);
  const TopTools_ListOfShape edge_collection = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
  const ShapeAnalysis_Edge edge_analyzer;
  for (auto edge = edge_collection.begin (); edge != edge_collection.end (); ++edge) {
    if (edge_analyzer.IsSeam (TopoDS::Edge (*edge), face))
      return true;
  }
  return false;
}

bool
t8_cad::t8_geom_edge_is_seam (const int edge_index, const int face_index) const
{
  const auto face = t8_cad::t8_geom_get_cad_face (face_index);
  const auto edge = t8_cad::t8_geom_get_cad_edge (edge_index);
  const ShapeAnalysis_Edge edge_analyzer;
  return edge_analyzer.IsSeam (edge, face);
}

void
t8_cad::t8_geom_get_parameter_of_vertex_on_edge (const int vertex_index, const int edge_index, double *edge_param,
                                                 std::optional<double> reference_edge_param) const
{
  T8_ASSERT (t8_cad::t8_geom_is_vertex_on_edge (vertex_index, edge_index));
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));

  /* If the edge is not closed or the user did not provide any reference coords we just query
  the parameter of the vertex on the edge. But if the edge is closed and the user provided
  said reference parameters it gets more sophisticated:*/

  const bool edge_is_closed = t8_cad::t8_geom_edge_is_closed (edge_index);
  if (reference_edge_param && edge_is_closed) {
    /* Edge is closed and the user provided reference parameters. We iterate over all vertices of the edge.
    Since the edge is closed, the start and end vertex should be in the same physical location.
    We choose the point where the parameters are closer to the reference parameters. For debugging reasons
    we also check, if the points are really in the same physical location. */

    bool first_point = true;
    for (TopExp_Explorer dora (edge, TopAbs_VERTEX); dora.More (); dora.Next ()) {
      const TopoDS_Vertex current_vertex = TopoDS::Vertex (dora.Current ());

#if T8_ENABLE_DEBUG
      /* Check of point is really the same. */
      const gp_Pnt debug_reference_point = BRep_Tool::Pnt (vertex);
      const gp_Pnt debug_current_point = BRep_Tool::Pnt (current_vertex);
      T8_ASSERT (debug_reference_point.Distance (debug_current_point) <= Precision::Confusion ());
#endif

      if (first_point) {
        *edge_param = BRep_Tool::Parameter (current_vertex, edge);
        first_point = false;
      }
      else {
        double other_param = BRep_Tool::Parameter (current_vertex, edge);
        if (std::abs (other_param - reference_edge_param.value ())
            < std::abs (*edge_param - reference_edge_param.value ()))
          *edge_param = other_param;
      }
    }
  }
  else {
    *edge_param = BRep_Tool::Parameter (vertex, edge);
  }
}

void
t8_cad::t8_geom_get_parameters_of_vertex_on_face (const int vertex_index, const int face_index, double face_params[2],
                                                  std::optional<std::span<const double, 2>> reference_face_params) const
{
  /* DISCLAIMER: This function is overly complicated and I do not understand why the simpler versions to not work.
  The overly complicated part is only the edge case where the vertex lies on a seam and reference parameters are provided.
  The (simpler) plan was as follows: We have the vertex and face. But in the topology of the face, the vertex exists two times,
  on the one side of the face and on the other. Bot sides are on the same physical location, because it is a seam. Thats
  why both vertices have the same coordinates. But they do not have the same parameters on the surface. The simple plan
  is to use a TopExp_Explorer to iterate over all vertices of the face and to look at all vertices with the right coords.
  From all those vertices we choose the one, which has the closest parameters to our reference parameters.
  But for whatever reason the surface in my test case only has vertices with one parameterset. Not with the other set.
  And thats why this algorithm does not work. What you see here is the more complicated workaround:
  */

  T8_ASSERT (t8_cad::t8_geom_is_vertex_on_face (vertex_index, face_index));
  std::optional<gp_Pnt2d> uv = std::nullopt; /** Final parameters on the surface */
  TopoDS_Vertex vertex = TopoDS::Vertex (cad_shape_vertex_map.FindKey (vertex_index));
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));

  /* If the surface is not closed or the user did not provide any reference params we just query the parameters. */
  if (!t8_cad::t8_geom_surface_is_closed (face_index) || !reference_face_params) {
    uv.emplace (BRep_Tool::Parameters (vertex, face));
  }

  /* If the vertex is not on the seam we can also just query the parameters. */
  else {
    const bool is_on_seam = t8_cad::t8_geom_vertex_is_on_seam (vertex_index, face_index);
    if (!is_on_seam) {
      uv.emplace (BRep_Tool::Parameters (vertex, face));
    }
    /* Now the overcomplicated edge case: The user provided reference params and the vertex is on the seam of the surface. */
    else {
      /* Convert reference parameters to a OCCT point */
      gp_Pnt2d reference_point (reference_face_params.value ()[0], reference_face_params.value ()[1]);

      /* The vertex is on a seam and so we have to get all vertices and check which one is closer.
      But since the vertices we are searching for are not in the topology of the face we have to take a detour.
      We iterate over all edges of the surface to find the seam edge. Yes, the seam EDGE. Not the seam EDGES.
      I would think that there are two seam edges. On both sides of the face there is one. And for a cylindrical
      surface the TopExp_Explorer gives us four edges. Two circles and two lines. I would say, that both lines are
      seams, because they are the same line. But OCCT only marks one of the two as seam. We now try to find this seam.
      */

      /** List of all edges connected to the vertex */
      const TopTools_ListOfShape vertex_edges = cad_shape_vertex2edge_map.FindFromIndex (vertex_index);
      const ShapeAnalysis_Edge edge_analyzer;
      std::optional<TopoDS_Edge> seam = std::nullopt;

      for (TopExp_Explorer dora (face, TopAbs_EDGE); dora.More (); dora.Next ()) {
        const TopoDS_Edge current_edge = TopoDS::Edge (dora.Current ());
        /* The edge has to be a seam and has to be connected to our vertex. */
        const bool current_edge_connected_to_vertex = vertex_edges.Contains (current_edge);
        const bool current_edge_is_seam = edge_analyzer.IsSeam (current_edge, face);
        if (current_edge_connected_to_vertex && current_edge_is_seam) {
          seam.emplace (current_edge);
          break;
        }
      }
      /* Hopefully we found a seam. If not something has gone catastrophically wrong or I overlooked something.
      Probably the latter. */
      SC_CHECK_ABORTF (seam.has_value (), "Error: Could not find seam of periodic surface.");

      /* We now have the seam and can iterate over all edges of the face again. Every edge matching our seam is also a seam (in my eyes).
      In the case of our cylinder this would be the two lines. */
      for (TopExp_Explorer dora (face, TopAbs_EDGE); dora.More (); dora.Next ()) {
        const TopoDS_Edge current_edge = TopoDS::Edge (dora.Current ());
        if (seam->IsSame (current_edge)) {
          /* If the edge matches our seam we can query the parameter of the vertex on the edge on the face.
          Why also on the face and not only on the edge? I really don't know but I had all information to call this function and it seems to work. */
          Standard_Real first, last;
          const double curve_param = BRep_Tool::Parameter (vertex, current_edge, face);
          const Handle_Geom2d_Curve pcurve_on_surface = BRep_Tool::CurveOnSurface (current_edge, face, first, last);
          /* We can now insert the queried parameter on the curve into the pcurve (the curve on the surface).
          If this is the first point we found uv is empty and we fill it. Otherwise, we define another point to
          stor our parameters and we save the parameters in iv which are closer to the reference parameters. */
          if (!uv) {
            uv.emplace ();
            pcurve_on_surface->D0 (curve_param, *uv);
          }
          /* If it is the second one we check which point is close. */
          else {
            gp_Pnt2d other_point;
            pcurve_on_surface->D0 (curve_param, other_point);
            if (other_point.Distance (reference_point) < uv->Distance (reference_point))
              uv = other_point;
          }
        }
      }
    }
  }

  /* Lastly, we put the parameters into our output array. */
  face_params[0] = uv->X ();
  face_params[1] = uv->Y ();
}

void
t8_cad::t8_geom_edge_parameter_to_face_parameters (
  const int edge_index, const int face_index, const double edge_param, double face_params_out[2],
  std::optional<std::span<const double, 2>> reference_face_params) const
{
  T8_ASSERT (t8_cad::t8_geom_is_edge_on_face (edge_index, face_index));
  Standard_Real first, last;
  gp_Pnt2d uv;
  TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (face_index));
  TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (edge_index));

  const bool is_seam = t8_cad::t8_geom_edge_is_seam (edge_index, face_index);
  if (is_seam && reference_face_params) {
    /* Convert reference parameters to OCCT point */
    gp_Pnt2d reference_point (reference_face_params.value ()[0], reference_face_params.value ()[1]);

    /* If the edge is a seam we have to get both points and check which one is closer.
       We iterate over all edges of the face and check if the edges are the same as the
       input edge. Then we check if the converted parameters are closer to the
       already computed parameters. */
    bool first_point = true;
    for (TopExp_Explorer dora (face, TopAbs_EDGE); dora.More (); dora.Next ()) {
      const TopoDS_Edge current_edge = TopoDS::Edge (dora.Current ());
      /* Check if edge is one of the seams. */
      if (edge.IsSame (current_edge)) {
        Handle_Geom2d_Curve curve_on_surface = BRep_Tool::CurveOnSurface (edge, face, first, last);
        /* If it is the first seam we compute the parameters. */
        if (first_point) {
          curve_on_surface->D0 (edge_param, uv);
          first_point = false;
        }
        /* If it is the second one we check which point is close. */
        else {
          gp_Pnt2d other_point;
          curve_on_surface->D0 (edge_param, other_point);
          if (other_point.Distance (reference_point) < uv.Distance (reference_point))
            uv = other_point;
        }
      }
    }
  }
  else {
    /* If the edge is not a seam or if no reference parameters were provided we just
       use the normal edge. */
    Handle_Geom2d_Curve curve_on_surface = BRep_Tool::CurveOnSurface (edge, face, first, last);
    curve_on_surface->D0 (edge_param, uv);
  }
  face_params_out[0] = uv.X ();
  face_params_out[1] = uv.Y ();
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

bool
t8_cad::t8_geom_edge_is_closed (int edge_index) const
{
  const Handle_Geom_Curve cad_edge = t8_geom_get_cad_curve (edge_index);
  return cad_edge->IsClosed ();
}

bool
t8_cad::t8_geom_surface_is_closed (int geometry_index) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (geometry_index);
  return cad_surface->IsUClosed () || cad_surface->IsVClosed ();
}

bool
t8_cad::t8_geom_surface_is_closed (int geometry_index, int direction) const
{
  const Handle_Geom_Surface cad_surface = t8_geom_get_cad_surface (geometry_index);
  switch (direction) {
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
