/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/*
  This test checks if the t8_cad_boundary struct works properly.
  It reads a .brep file and its corresponding .msh file. After constructing the struct, the distance
  between the mesh node and the geometry is calculated. If the distance does not need exceed a certain tolerance
  it passes the test.
*/

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_boundary_node_list.hxx>
#include <t8_cmesh/t8_cmesh_cad/t8_cmesh_cad_boundary.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>
#include "t8_test_data_dir.h"

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp.hxx>
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
    std::string testfile_brep = std::string (T8_TEST_DATA_DIR) + brep_files[GetParam ()];
    cad_shape = read_brep_file (testfile_brep);
    TopExp::MapShapes (cad_shape, TopAbs_VERTEX, cad_shape_vertex_map);
    TopExp::MapShapes (cad_shape, TopAbs_EDGE, cad_shape_edge_map);
    TopExp::MapShapes (cad_shape, TopAbs_FACE, cad_shape_face_map);

    std::string testfile_msh = std::string (T8_TEST_DATA_DIR) + mesh_files[GetParam ()];
    cmesh = t8_cmesh_from_msh_file (testfile_msh.c_str (), 0, sc_MPI_COMM_WORLD, 3, 0, 0);
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

#define NUMBER_OF_TESTFILES 2

  static constexpr std::array<const char*, NUMBER_OF_TESTFILES> mesh_files = { "/simple_test_case", "/D150" };

  static constexpr std::array<const char*, NUMBER_OF_TESTFILES> brep_files = { "/simple_test_case", "/D150" };

  int file;
};

TEST_P (t8_gtest_cad_boundary, geom_data_map_test)
{
  constexpr double tolerance = 1e-6;

  t8_boundary_node_geom_data_map boundary_node_map = t8_boundary_node_geom_data_map (cad_shape, cmesh, tolerance);
  const std::unordered_map<t8_gloidx_t, t8_geom_data>& geom_data_map
    = boundary_node_map.get_boundary_node_geom_data_map ();
  t8_productionf ("geom_data_map created with the size %lu\n", geom_data_map.size ());
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
      /* Check that entity_tag is a valid index in the edge map */
      ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
      ASSERT_LE (iter.second.entity_tag, cad_shape_vertex_map.Extent ()) << "entity_tag index too high";

      const gp_Pnt point = BRep_Tool::Pnt (TopoDS::Vertex (cad_shape_vertex_map.FindKey (iter.second.entity_tag)));

      const double dist = point.Distance (gp_Pnt (cmesh_x_coords, cmesh_y_coords, cmesh_z_coords));

      /* Check that distance between mesh node and geometry vertex is less or equal than tolerance */
      EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;
    }

    /* For curves of geometry */
    else if (iter.second.entity_dim == 1) {
      Standard_Real first, last;

      /* Check that entity_tag is a valid index in the edge map */
      ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
      ASSERT_LE (iter.second.entity_tag, cad_shape_edge_map.Extent ()) << "entity_tag index too high";

      TopoDS_Edge edge = TopoDS::Edge (cad_shape_edge_map.FindKey (iter.second.entity_tag));
      Handle (Geom_Curve) curve = BRep_Tool::Curve (edge, first, last);

      ASSERT_FALSE (curve.IsNull ()) << "Null curve handle for edge " << iter.second.entity_tag;

      /* Reverse projection */
      curve->D0 (iter.second.location_on_curve[0], point);

      const double dist = point.Distance (gp_Pnt (cmesh_x_coords, cmesh_y_coords, cmesh_z_coords));

      /* Check that distance between mesh node and geometry curve is less or equal than tolerance */
      EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;
    }

    else if (iter.second.entity_dim == 2) {
      /* Check that entity_tag is a valid index in the edge map */
      ASSERT_GE (iter.second.entity_tag, 1) << "entity_tag index too low";
      ASSERT_LE (iter.second.entity_tag, cad_shape_face_map.Extent ()) << "entity_tag index too high";

      TopoDS_Face face = TopoDS::Face (cad_shape_face_map.FindKey (iter.second.entity_tag));
      Handle_Geom_Surface surface = BRep_Tool::Surface (face);

      ASSERT_FALSE (surface.IsNull ()) << "Null curve handle for edge " << iter.second.entity_tag;

      /* Reverse projection */
      surface->D0 (iter.second.location_on_curve[0], iter.second.location_on_curve[1], point);

      const double dist = point.Distance (gp_Pnt (cmesh_x_coords, cmesh_y_coords, cmesh_z_coords));

      /* Check that distance between mesh node and geometry curve is less or equal than tolerance */
      EXPECT_LE (dist, tolerance) << "Distance exceeds tolerance for edge " << iter.second.entity_tag;
    }
  }
}
INSTANTIATE_TEST_SUITE_P (t8_gtest_geom_data_map, t8_gtest_cad_boundary, testing::Range (0, NUMBER_OF_TESTFILES));
