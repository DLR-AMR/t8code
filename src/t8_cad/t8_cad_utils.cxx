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

#include <t8_cad/t8_cad_utils.hxx>

#if T8_WITH_OCC
#include <TopoDS.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepTools.hxx>

/* *INDENT-OFF* */

TopoDS_Shape
t8_cad_read_cad_file (const char * filename)
{
  std::string format;
  const std::string name = filename;
  const int idot = name.find_last_of(".");
  if (idot != (int)std::string::npos && idot < (int)(name.length() - 1)) {
    format = name.substr(idot);
  }
  else {
    SC_ABORTF ("Unable to parse CAD file format of %s", filename);
  }

  if (format == ".brep" || format == ".BREP") {
    t8_global_productionf ("Reading in BREP file %s \n", filename);
    TopoDS_Shape  shape;
    BRep_Builder  builder;
    std::ifstream is (name);
    BRepTools::Read (shape, is, builder);
    is.close ();
    if (shape.IsNull ()) {
      SC_ABORTF ("Could not read BREP file or BREP file contains no shape.");
    }
    return shape;
  }
  else if (format == ".step" || format == ".STEP"
           || format == ".stp" || format == ".STP") {
    t8_global_productionf ("Reading in STEP file %s \n", filename);
    STEPControl_Reader reader;
    if (reader.ReadFile(filename) != IFSelect_RetDone) {
      SC_ABORTF ("Could not read STEP file %s", filename);
    }
#if T8_ENABLE_DEBUG
    reader.PrintCheckLoad(0, IFSelect_ItemsByEntity);
#else
    reader.PrintCheckLoad(1, IFSelect_ItemsByEntity);
#endif
    reader.NbRootsForTransfer();
    reader.TransferRoots();
    return reader.OneShape();
  }
  else if (format == ".iges" || format == ".IGES"
           || format == ".igs" || format == ".IGS") {
    t8_global_productionf ("Reading in IGES file %s \n", filename);
    IGESControl_Reader reader;
    if (reader.ReadFile(filename) != IFSelect_RetDone) {
      SC_ABORTF ("Could not read IGES file %s", filename);
    }
#if T8_ENABLE_DEBUG
    reader.PrintCheckLoad(0, IFSelect_ItemsByEntity);
#else
    reader.PrintCheckLoad(1, IFSelect_ItemsByEntity);
#endif
    reader.NbRootsForTransfer();
    reader.TransferRoots();
    return reader.OneShape();
  }
  else {
    SC_ABORTF ("Error: Unknown CAD file format: %s", format.c_str());
  }
}

TopoDS_Shape
t8_cad_make_axis_aligned_hex_element_shape (const double *vertex1,
                                            const double *vertex2)
{
  const gp_Pnt pnt1{vertex1[0], vertex1[1], vertex1[2]};
  const gp_Pnt pnt2{vertex2[0], vertex2[1], vertex2[2]};
  return BRepPrimAPI_MakeBox(pnt1, pnt2);
}

TopoDS_Shape
t8_cad_make_element_shape (const double *vertices, const t8_eclass_t eclass)
{
  TopoDS_Shape shape;
  const int num_vertices = t8_eclass_num_vertices[eclass];
  /* Create a point of each element corner */
  gp_Pnt *pnts = T8_ALLOC (gp_Pnt, num_vertices);
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex)
  {
    pnts[i_vertex] = gp_Pnt{
      vertices[i_vertex * 3],
      vertices[i_vertex * 3 + 1],
      vertices[i_vertex * 3 + 2]
    };
  }
  switch (t8_eclass_to_dimension[eclass])
  {
  case 0:
  {
    shape = BRepBuilderAPI_MakeVertex(pnts[0]);
    break;
  }
  case 1:
  {
    shape = BRepBuilderAPI_MakeEdge(pnts[0], pnts[1]);
    break;
  }
  case 2:
  {
    /* We create a wire out of each edge of the 2D element. */
    const int num_edges = t8_eclass_num_edges[eclass];
    /* Create a list fot the edges */
    TopTools_ListOfShape edge_list;
    for (int i_edge = 0; i_edge < num_edges; ++i_edge)
    {
      /* Get edge vertices */
      const int vertex1 = t8_face_vertex_to_tree_vertex[eclass][0][0];
      const int vertex2 = t8_face_vertex_to_tree_vertex[eclass][0][1];
      /* Make an edge shape from the vertices and append it to the list */
      edge_list.Append (BRepBuilderAPI_MakeEdge (pnts[vertex1],
                                                 pnts[vertex2]));
    }
    /* There should be num_edges edges in the list */
    T8_ASSERT (edge_list.Size() == num_edges);
    /* Make a wire and then a face shape out of the wire */
    BRepBuilderAPI_MakeWire mkwire{};
    mkwire.Add (edge_list);
    shape = BRepBuilderAPI_MakeFace (mkwire.Wire(), true);
    break;
  }
  case 3:
  {
    /* Make a face shape for each element face and sew them together */
    const int num_faces = t8_eclass_num_faces[eclass];
    BRepBuilderAPI_Sewing sewer{Precision::Confusion(), true, false, false, false};
    for (int i_face = 0; i_face < num_faces; ++i_face)
    {
      /* Create a wire for each face */
      TopoDS_Wire wire {};
      const t8_eclass_t face_eclass = (t8_eclass_t)t8_eclass_face_types[eclass][i_face];
      const int num_face_edges = t8_eclass_num_edges[face_eclass];
      /* Create a list fot the edges */
      TopTools_ListOfShape edge_list;
      for (int i_face_edge = 0; i_face_edge < num_face_edges; ++i_face_edge)
      {
        /* Get the vertices of the edge of the current face */
        const int face_vertex1 = t8_face_vertex_to_tree_vertex[face_eclass][i_face_edge][0];
        const int face_vertex2 = t8_face_vertex_to_tree_vertex[face_eclass][i_face_edge][1];
        const int tree_vertex1 = t8_face_vertex_to_tree_vertex[eclass][i_face][face_vertex1];
        const int tree_vertex2 = t8_face_vertex_to_tree_vertex[eclass][i_face][face_vertex2];
        /* Make an edge shape out of the vertices and add it to the list */
        edge_list.Append (BRepBuilderAPI_MakeEdge (pnts[tree_vertex1],
                                                   pnts[tree_vertex2]));
      }
      /* There should be num_face_edges edges in the list */
      T8_ASSERT (edge_list.Size() == num_face_edges);
      /* Make a wire and then a face shape out of the wire and add it to the sewer */
      BRepBuilderAPI_MakeWire mkwire{};
      mkwire.Add (edge_list);
      const TopoDS_Face face = BRepBuilderAPI_MakeFace (mkwire.Wire(), true);

      sewer.Add (face);
    }
    sewer.Perform();
    /* There should remain no unsewed edges and the shape should be a shell */
    if (sewer.NbFreeEdges())
    {
      BRepTools::Write(sewer.SewedShape(), "fail.brep");
    }
    T8_ASSERT (sewer.NbFreeEdges() == 0);
    T8_ASSERT (sewer.SewedShape().ShapeType() == TopAbs_SHELL);
    shape = BRepBuilderAPI_MakeSolid (TopoDS::Shell(sewer.SewedShape()));
    break;
  }
  default:
    break;
  }
  T8_FREE (pnts);
  return shape;
}

/* *INDENT-OFF* */

#endif /* T8_WITH_OCC */
