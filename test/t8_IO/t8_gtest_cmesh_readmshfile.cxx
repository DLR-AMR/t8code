/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <gtest/gtest.h>
#include <unistd.h> /* Needed to check for file access */
#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_trees.h>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>
#include "t8_test_data_dir.h"
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <algorithm>

/* In this file we test the msh file (gmsh) reader of the cmesh.
 * Currently, we support version 2 and 4 ascii.
 * We read a mesh file and check whether the constructed cmesh is correct.
 * We also try to read version 2 binary and version 4 binary
 * formats. All are not supported and we expect the reader to catch this.
 */

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_ascii)
{

  auto fileprefix = t8_test_data_dir / "test_msh_file_vers4_ascii";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix.c_str ());

  t8_debugf ("Checking msh file version 4 ascii...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_from_msh_file (&cmesh, fileprefix.c_str (), 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from ascii version 4, but should be able to.";

  /* The cmesh was read successfully and we need to destroy it. */
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_readmshfile, test_msh_file_vers4_bin)
{

  auto fileprefix = t8_test_data_dir / "test_msh_file_vers4_bin";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix.c_str ());

  t8_debugf ("Checking msh file version 4 binary...\n");

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_from_msh_file (&cmesh, fileprefix.c_str (), 1, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh == NULL) << "Expected fail of reading binary msh file v.4, but did not fail.";

  t8_debugf ("Error handling successful.\n");
}

/* Reading a msh file with a dimension it has no elements of must fail
 * cleanly while reading the elements and return a NULL cmesh. */
TEST (t8_cmesh_readmshfile, test_msh_file_vers4_ascii_wrong_dim)
{
  auto fileprefix = t8_test_data_dir / "test_msh_file_vers4_ascii";
  char filename[BUFSIZ];

  snprintf (filename, BUFSIZ, "%s.msh", fileprefix.c_str ());

  ASSERT_FALSE (access (filename, R_OK)) << "Could not open file " << filename;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  /* The file only contains 2D elements. */
  t8_cmesh_from_msh_file (&cmesh, fileprefix.c_str (), 1, sc_MPI_COMM_WORLD, 3, 0, 0);
  ASSERT_TRUE (cmesh == NULL) << "Expected fail of reading 3D elements from a 2D msh file, but did not fail.";
}

/* Returns the boundary condition that the test files with physical groups assign to the face with the
 * given vertex coordinates. The files mesh the unit square and the unit cube. Their boundary at x = 0 is
 * labeled "inflow", the one at x = 1 "outflow" and all other boundaries \a wall_label. Faces in the
 * interior do not belong to a physical group and thus have an empty label. */
static std::string_view
t8_test_msh_expected_boundary_condition (const std::vector<std::array<double, 3>> &face_coordinates, const int dim,
                                         const std::string_view wall_label)
{
  /* Checks whether all vertices of the face lie on the plane where the coordinate \a axis equals \a value.
   * Gmsh writes the coordinates of nodes on the boundary of the unit square/cube exactly as 0 or 1,
   * so we can compare without a tolerance. */
  const auto all_on_plane = [&face_coordinates] (const int axis, const double value) {
    return std::all_of (
      face_coordinates.begin (), face_coordinates.end (),
      [axis, value] (const std::array<double, 3> &coordinates) { return coordinates[axis] == value; });
  };
  /* The boundaries in x-direction have their own physical groups. */
  if (all_on_plane (0, 0)) {
    return "inflow";
  }
  if (all_on_plane (0, 1)) {
    return "outflow";
  }
  /* All other boundaries belong to the wall group. In 2D these are the boundaries in y-direction,
   * in 3D the boundaries in y- and z-direction. */
  for (int axis = 1; axis < dim; ++axis) {
    if (all_on_plane (axis, 0) || all_on_plane (axis, 1)) {
      return wall_label;
    }
  }
  /* The face is not on the boundary of the domain. */
  return "";
}

/* Read a msh file with physical groups on the boundary and check that each tree face
 * got the name of its physical group as boundary condition.
 * The expected boundary condition of a face is derived from the position of its vertices,
 * see t8_test_msh_expected_boundary_condition. */
static void
t8_test_msh_file_boundary_conditions (const char *filename, const int dim, const std::string_view wall_label)
{
  /* Read the msh file into a replicated cmesh. */
  const auto fileprefix = t8_test_data_dir / filename;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_from_msh_file (&cmesh, fileprefix.c_str (), 0, sc_MPI_COMM_WORLD, dim, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from " << fileprefix;

  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t ltree_id = 0; ltree_id < num_local_trees; ++ltree_id) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, ltree_id);
    /** The coordinates of the tree vertices, three values per vertex. */
    const double *tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltree_id);
    /* Every tree must have exactly one boundary condition per face, including its interior faces. */
    const t8_boundary_conditions<std::string_view> boundary_conditions
      = t8_cmesh_get_boundary_conditions (cmesh, ltree_id);
    ASSERT_EQ (boundary_conditions.size (), static_cast<size_t> (t8_eclass_num_faces[eclass]));
    for (int face = 0; face < t8_eclass_num_faces[eclass]; ++face) {
      /* Collect the coordinates of the face vertices to determine where the face lies. */
      const int num_face_vertices = t8_eclass_num_vertices[t8_eclass_face_types[eclass][face]];
      std::vector<std::array<double, 3>> face_coordinates (num_face_vertices);
      for (int i_vertex = 0; i_vertex < num_face_vertices; ++i_vertex) {
        const int tree_vertex = t8_face_vertex_to_tree_vertex[eclass][face][i_vertex];
        std::copy_n (tree_vertices + 3 * tree_vertex, 3, face_coordinates[i_vertex].begin ());
      }
      /* Compare the boundary condition read from the file with the one expected at this position. */
      EXPECT_EQ (boundary_conditions[face], t8_test_msh_expected_boundary_condition (face_coordinates, dim, wall_label))
        << "Wrong boundary condition at face " << face << " of tree " << ltree_id << ".";
    }
  }
  t8_cmesh_destroy (&cmesh);
}

/* The 2D file meshes the unit square with triangles. The top and bottom boundary
 * belong to the physical group 7, which has no name. */
TEST (t8_cmesh_readmshfile, test_msh_file_boundary_conditions_2d)
{
  t8_test_msh_file_boundary_conditions ("test_msh_file_boundary_conditions_2d", 2, "physical_group_7");
}

/* The 3D file meshes the unit cube with hexahedra. All boundaries except
 * x = 0 and x = 1 belong to the physical group "wall". */
TEST (t8_cmesh_readmshfile, test_msh_file_boundary_conditions_3d)
{
  t8_test_msh_file_boundary_conditions ("test_msh_file_boundary_conditions_3d", 3, "wall");
}

/* A file without physical groups must create a cmesh without boundary conditions. */
TEST (t8_cmesh_readmshfile, test_msh_file_without_boundary_conditions)
{
  /* This file has an $Entities section, but none of its entities belongs to a physical group. */
  const auto fileprefix = t8_test_data_dir / "test_msh_file_vers4_ascii";
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_from_msh_file (&cmesh, fileprefix.c_str (), 0, sc_MPI_COMM_WORLD, 2, 0, 0);
  ASSERT_TRUE (cmesh != NULL) << "Could not read cmesh from " << fileprefix;
  /* The boundary condition getters abort on a cmesh without boundary conditions.
   * Hence, we check that no boundary condition handler was created instead. */
  EXPECT_EQ (detail::t8_cmesh_get_boundary_condition_handler (cmesh), nullptr);
  t8_cmesh_destroy (&cmesh);
}
