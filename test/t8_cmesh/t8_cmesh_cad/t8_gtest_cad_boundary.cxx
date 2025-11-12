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
#include <t8_cmesh.h>
#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_trees.h>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_cad_boundary.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp.hxx>
#include <t8_cmesh_readmshfile.h>
#include <TopTools_IndexedMapOfShape.hxx>

TopoDS_Shape
read_brep_file (std::string fileprefix)
{
  TopoDS_Shape cad_shape;
  BRep_Builder bob;
  std::string current_file (fileprefix);
  std::ifstream is (current_file + ".brep");
  if (is.is_open () == false) {
    SC_ABORTF ("Cannot find the file %s.brep.\n", fileprefix.c_str ());
  }
  BRepTools::Read (cad_shape, is, bob);
  is.close ();
  if (cad_shape.IsNull ()) {
    SC_ABORTF ("Could not read brep file or brep file contains no shape.");
  }

  return cad_shape;
}

class t8_gtest_cad_boundary: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
    cad_shape = read_brep_file (brep_files[GetParam ()]);
    TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
    TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
    TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);

    cmesh = t8_cmesh_from_msh_file (mesh_files[GetParam ()], 0, sc_MPI_COMM_WORLD, 3, 0, 0);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  TopoDS_Shape cad_shape;
  TopTools_IndexedMapOfShape cad_shape_vertex_map; /**< Map of all TopoDS_Vertex in shape. */
  TopTools_IndexedMapOfShape cad_shape_edge_map;   /**< Map of all TopoDS_Edge in shape. */
  TopTools_IndexedMapOfShape cad_shape_face_map;   /**< Map of all TopoDS_Face in shape. */
  t8_cmesh_t cmesh;

  const char* mesh_files[7] /* .msh file prefixes */
    = { "test/testfiles/simple_test_case_168",
        "test/testfiles/simple_test_case_288",
        "test/testfiles/simple_test_case_1414",
        "test/testfiles/simple_test_case_5147",
        "test/testfiles/simple_test_case_30596",
        "test/testfiles/D150_1791",
        "test/testfiles/D150_5129" };

  const char* brep_files[7] /* .brep file prefixes */
    = { "test/testfiles/simple_test_case",
        "test/testfiles/simple_test_case",
        "test/testfiles/simple_test_case",
        "test/testfiles/simple_test_case",
        "test/testfiles/simple_test_case",
        "test/testfiles/D150",
        "test/testfiles/D150" };

  int file;
};

TEST_P (t8_gtest_cad_boundary, geom_data_map_test)
{
  double tolerance = 1e-6;

  t8_boundary_node_geom_data_map boundary_node_map = t8_boundary_node_geom_data_map (cad_shape, cmesh, tolerance);
  std::unordered_map geom_data_map = boundary_node_map.get_boundary_node_geom_data_map ();
  t8_productionf ("Geom Data Map Created with the size %lu\n", geom_data_map.size ());
  ASSERT_EQ (cmesh->boundary_node_list->get_boundary_node_list ().size (), geom_data_map.size ());

  for (auto iter : geom_data_map) {
    gp_Pnt point;

    /* Get coordinates of mesh node */
    const tree_vertex_list tree_list = cmesh->vertex_connectivity->vertex_to_trees (iter.first);
    t8_locidx_t local_tree_id = tree_list.at (0).first;
    int local_vertex_id = tree_list.at (0).second;
    double* vertices = (double*) t8_cmesh_get_tree_vertices (cmesh, local_tree_id);
    const double cmesh_x_coords = vertices[3 * local_vertex_id];
    const double cmesh_y_coords = vertices[3 * local_vertex_id + 1];
    const double cmesh_z_coords = vertices[3 * local_vertex_id + 2];

    /* For vertices of geometry */
    if (iter.second.entity_dim == 0) {
      try {
        /* Check that entity_tag is a valid index in the edge map */
        ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
        ASSERT_LE (iter.second.entity_tag, cad_shape_vertex_map.Extent ()) << "entity_tag index too high";

        point = BRep_Tool::Pnt (TopoDS::Vertex (cad_shape_vertex_map.FindKey (iter.second.entity_tag)));

        double dx = cmesh_x_coords - point.X ();
        double dy = cmesh_y_coords - point.Y ();
        double dz = cmesh_z_coords - point.Z ();

        double dist = std::sqrt (dx * dx + dy * dy + dz * dz);

        /* Check that distance between mesh node and geometry vertex is less or equal than tolerance */
        EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;

      } catch (const Standard_Failure& e) { /* OCC Exception Handling*/
        FAIL () << "Open CASCADE exception: " << e.GetMessageString ();
      } catch (...) {
        FAIL () << "Unknown exception caught during curve evaluation";
      }
    }

    /* For curves of geometry */
    else if (iter.second.entity_dim == 1) {
      try {
        Standard_Real first, last;

        /* Check that entity_tag is a valid index in the edge map */
        ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
        ASSERT_LE (iter.second.entity_tag, cad_shape_edge_map.Extent ()) << "entity_tag index too high";

        TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (iter.second.entity_tag));
        Handle (Geom_Curve) curve = BRep_Tool::Curve (edge, first, last);

        ASSERT_FALSE (curve.IsNull ()) << "Null curve handle for edge " << iter.second.entity_tag;

        /* Reverse projection */
        curve->D0 (iter.second.location_on_curve[0], point);

        double dx = cmesh_x_coords - point.X ();
        double dy = cmesh_y_coords - point.Y ();
        double dz = cmesh_z_coords - point.Z ();

        double dist = std::sqrt (dx * dx + dy * dy + dz * dz);

        /* Check that distance between mesh node and geometry curve is less or equal than tolerance */
        EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;

      } catch (const Standard_Failure& e) { /* OCC Exception Handling*/
        FAIL () << "Open CASCADE exception: " << e.GetMessageString ();
      } catch (...) {
        FAIL () << "Unknown exception caught during curve evaluation";
      }
    }

    else if (iter.second.entity_dim == 2) {
      try {
        /* Check that entity_tag is a valid index in the edge map */
        ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
        ASSERT_LE (iter.second.entity_tag, cad_shape_face_map.Extent ()) << "entity_tag index too high";

        TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (iter.second.entity_tag));
        Handle_Geom_Surface surface = BRep_Tool::Surface (face);

        ASSERT_FALSE (surface.IsNull ()) << "Null curve handle for edge " << iter.second.entity_tag;

        /* Reverse projection */
        surface->D0 (iter.second.location_on_curve[0], iter.second.location_on_curve[1], point);

        double dx = cmesh_x_coords - point.X ();
        double dy = cmesh_y_coords - point.Y ();
        double dz = cmesh_z_coords - point.Z ();

        double dist = std::sqrt (dx * dx + dy * dy + dz * dz);

        /* Check that distance between mesh node and geometry curve is less or equal than tolerance */
        EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;

      } catch (const Standard_Failure& e) { /* OCC Exception Handling*/
        FAIL () << "Open CASCADE exception: " << e.GetMessageString ();
      } catch (...) {
        FAIL () << "Unknown exception caught during curve evaluation";
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_geom_data_map, t8_gtest_cad_boundary, testing::Range (0, 7));
