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

#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>
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
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Shape.hxx>
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

class t8_gtest_cad_boundary: public testing::Test {
 protected:
  void
  SetUp () override
  {
    cad_shape = read_brep_file ("test/testfiles/coladose");
    cmesh = t8_cmesh_from_msh_file ("test/testfiles/coladose", 0, sc_MPI_COMM_WORLD, 3, 0, 0);
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
  }

  TopoDS_Shape cad_shape;
  t8_cmesh_t cmesh;
};

TEST_F (t8_gtest_cad_boundary, some_random_ass_name)
{
  t8_boundary_node_geom_data_map boundary_node_map = t8_boundary_node_geom_data_map (cad_shape, cmesh, 1e-6);
  std::unordered_map geom_data_map = boundary_node_map.get_boundary_node_geom_data_map ();
  t8_debugf ("Geom Data Map Created with the size %lu\n", geom_data_map.size ());
  for (auto iter : geom_data_map) {
    if (iter.second.entity_dim == 2) {
      TopoDS_Face surface = BRep_Tool::Surface (TopoDS::Face (cad_shape_face_map.FindKey (iter.second.entity_tag)));
    }
  }
}
