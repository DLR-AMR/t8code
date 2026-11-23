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

/** \file t8_gtest_cmesh_cgns_reader.cxx
 * Tests for the CGNS reader. The files in the test data directory are read and the resulting cmeshes
 * are checked.
 *
 * The files test_cgns_*.cgns were written with the CGNS library (version 4.5.2, HDF5 format):
 *  - test_cgns_two_zones.cgns: Two unstructured zones of hexahedra.
 *      Zone A: two hexahedra [0,1]x[0,1]x[0,1] and [1,2]x[0,1]x[0,1] (section HEXA_8),
 *              boundary condition "wall" (BCWall, family "WallFamily", GridLocation FaceCenter on a QUAD_4
 *              section) on the faces at y = 0, Vertex connectivity "AtoB" to zone B at x = 2.
 *      Zone B: one hexahedron [2,3]x[0,1]x[0,1], boundary condition "top" (BCOutflow, GridLocation Vertex)
 *              on the face at z = 1, periodic FaceCenter connectivity "BtoA" from x = 3 to x = 0 of zone A
 *              with Translation (-3, 0, 0).
 *  - test_cgns_two_structured_zones.cgns: The same mesh with structured zones.
 *      Zone A: 3 x 2 x 2 vertices, boundary condition "wall" as PointRange on the plane j = 1,
 *              GridConnectivity1to1 "AtoB" with Transform [3, 1, 2].
 *      Zone B: 2 x 2 x 2 vertices with the index axes cyclically permuted: i runs along y, j along z and k
 *              along x. Periodic GridConnectivity1to1 "BtoA" from k = 2 to i = 1 of zone A with
 *              Transform [2, 3, 1], boundary condition "top" as PointRange on the plane j = 2.
 *  - test_cgns_mixed_2d.cgns: One unstructured 2D zone with a MIXED section of the quadrilateral
 *      (0,0),(1,0),(1,1),(0,1) and the triangle (2,0.5),(1,1),(1,0) and a BAR_2 section with the edge at y = 0,
 *      which carries the boundary condition "bottom" (BCWall, GridLocation EdgeCenter, PointRange).
 *
 * The other files are meshes of the DLR flow solver TRACE:
 *  - simpleElements_1-1.cgns: One unstructured zone of two hexahedra with periodic FaceCenter
 *      connectivities in two directions (the domain is one cell wide) and two wall boundary conditions.
 *  - box_9x9x9_peroidicY.cgns: One structured zone of 9 x 9 x 9 hexahedra with periodic
 *      GridConnectivity1to1 nodes in the j direction and four boundary conditions as PointRange.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_cgns_reader.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/t8_cmesh_boundary_conditions.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include "t8_test_data_dir.h"

#include <cgnslib.h>

#include <set>
#include <string>
#include <unistd.h> /* Needed to check for file access */

/** Read a file of the test data directory with the reader.
 * \return The cmesh. The test is aborted if the file is missing or cannot be read.
 */
static t8_cmesh_t
t8_cgns_read_test_file (const std::string &name, t8_cmesh_cgns_reader &reader)
{
  const std::string filename = (t8_test_data_dir / name).string ();
  EXPECT_FALSE (access (filename.c_str (), R_OK)) << "Could not open file " << filename;
  reader = t8_cmesh_cgns_reader (filename);
  t8_cmesh_t cmesh = reader.read (sc_MPI_COMM_WORLD);
  EXPECT_NE (cmesh, nullptr) << "Could not read " << filename;
  return cmesh;
}

/** Check that \a ltree is connected at \a face to \a expected_neighbor at \a expected_dual_face. */
static void
t8_cgns_expect_neighbor (t8_cmesh_t cmesh, const t8_locidx_t ltree, const int face, const t8_locidx_t expected_neighbor,
                         const int expected_dual_face, const int expected_orientation)
{
  EXPECT_FALSE (t8_cmesh_tree_face_is_boundary (cmesh, ltree, face)) << "tree " << ltree << " face " << face;
  int dual_face = -1, orientation = -1;
  const t8_locidx_t neighbor = t8_cmesh_get_face_neighbor (cmesh, ltree, face, &dual_face, &orientation);
  EXPECT_EQ (neighbor, expected_neighbor) << "tree " << ltree << " face " << face;
  EXPECT_EQ (dual_face, expected_dual_face) << "tree " << ltree << " face " << face;
  EXPECT_EQ (orientation, expected_orientation) << "tree " << ltree << " face " << face;
}

/** Check the cmesh read from the two zone files.
 * \param [in] rotated_zone_b  True for the structured file, in which the index axes of zone B are permuted.
 */
static void
t8_cgns_check_two_zone_cmesh (t8_cmesh_t cmesh, const t8_cmesh_cgns_reader &reader, const bool rotated_zone_b)
{
  /* Trees and zones */
  EXPECT_EQ (reader.get_dimension (), 3);
  EXPECT_EQ (t8_cmesh_get_dimension (cmesh), 3);
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), 3);
  ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), 3);
  for (t8_locidx_t itree = 0; itree < 3; ++itree) {
    EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, itree), T8_ECLASS_HEX);
  }
  ASSERT_EQ (reader.get_zones ().size (), 2u);
  EXPECT_EQ (reader.get_zones ()[0].name, "ZoneA");
  EXPECT_EQ (reader.get_zones ()[0].first_tree, 0);
  EXPECT_EQ (reader.get_zones ()[0].num_trees, 2);
  EXPECT_EQ (reader.get_zones ()[0].num_vertices, 12);
  EXPECT_EQ (reader.get_zones ()[1].name, "ZoneB");
  EXPECT_EQ (reader.get_zones ()[1].first_tree, 2);
  EXPECT_EQ (reader.get_zones ()[1].num_trees, 1);
  EXPECT_EQ (t8_cmesh_get_cgns_zone (cmesh, 0), 0);
  EXPECT_EQ (t8_cmesh_get_cgns_zone (cmesh, 1), 0);
  EXPECT_EQ (t8_cmesh_get_cgns_zone (cmesh, 2), 1);

  /* Vertex ordering: the t8code hex vertex i has the coordinates of the i-th corner in z-order.
   * In the rotated zone B the x axis of the tree is the z-order's third axis and the z axis is its second. */
  for (t8_locidx_t itree = 0; itree < 3; ++itree) {
    const double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
    ASSERT_NE (vertices, nullptr);
    const bool rotated = rotated_zone_b && itree == 2;
    for (int ivertex = 0; ivertex < 8; ++ivertex) {
      const int dx = rotated ? (ivertex >> 2) & 1 : ivertex & 1;
      const int dy = rotated ? ivertex & 1 : (ivertex >> 1) & 1;
      const int dz = rotated ? (ivertex >> 1) & 1 : (ivertex >> 2) & 1;
      EXPECT_DOUBLE_EQ (vertices[3 * ivertex], itree + dx) << "tree " << itree << " vertex " << ivertex;
      EXPECT_DOUBLE_EQ (vertices[3 * ivertex + 1], dy) << "tree " << itree << " vertex " << ivertex;
      EXPECT_DOUBLE_EQ (vertices[3 * ivertex + 2], dz) << "tree " << itree << " vertex " << ivertex;
    }
  }

  /* Neighbors: within zone A, across the zone interface and across the periodic interface.
   * In the rotated zone B the x- face of the tree is face 4 and the x+ face is face 5, the z+ face is face 3. */
  const int b_face_xminus = rotated_zone_b ? 4 : 0;
  const int b_face_xplus = rotated_zone_b ? 5 : 1;
  const int b_face_zplus = rotated_zone_b ? 3 : 5;
  t8_cgns_expect_neighbor (cmesh, 0, 1, 1, 0, 0);
  t8_cgns_expect_neighbor (cmesh, 1, 0, 0, 1, 0);
  t8_cgns_expect_neighbor (cmesh, 1, 1, 2, b_face_xminus, 0);
  t8_cgns_expect_neighbor (cmesh, 2, b_face_xminus, 1, 1, 0);
  t8_cgns_expect_neighbor (cmesh, 2, b_face_xplus, 0, 0, 0);
  t8_cgns_expect_neighbor (cmesh, 0, 0, 2, b_face_xplus, 0);
  for (t8_locidx_t itree = 0; itree < 3; ++itree) {
    for (int face = 0; face < 6; ++face) {
      const bool is_x_face = itree < 2 ? face < 2 : (face == b_face_xminus || face == b_face_xplus);
      if (!is_x_face) {
        EXPECT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, itree, face)) << "tree " << itree << " face " << face;
      }
    }
  }

  /* Boundary conditions */
  ASSERT_EQ (reader.get_boundary_conditions ().size (), 2u);
  const t8_cmesh_cgns_boundary_condition &wall = reader.get_boundary_conditions ()[0];
  EXPECT_EQ (wall.name, "wall");
  EXPECT_EQ (wall.family, "WallFamily");
  EXPECT_EQ (wall.type, static_cast<int> (CGNS_ENUMV (BCWall)));
  EXPECT_EQ (wall.type_name, "BCWall");
  EXPECT_EQ (wall.zone, 0);
  EXPECT_EQ (wall.num_faces, 2);
  const t8_cmesh_cgns_boundary_condition &top = reader.get_boundary_conditions ()[1];
  EXPECT_EQ (top.name, "top");
  EXPECT_EQ (top.family, "");
  EXPECT_EQ (top.zone, 1);
  EXPECT_EQ (top.num_faces, 1);
  for (t8_locidx_t itree = 0; itree < 3; ++itree) {
    for (int face = 0; face < 6; ++face) {
      std::string expected;
      if (itree < 2 && face == 2) {
        expected = "wall";
      }
      if (itree == 2 && face == b_face_zplus) {
        expected = "top";
      }
      EXPECT_EQ (t8_cmesh_get_boundary_condition (cmesh, itree, face), expected)
        << "tree " << itree << " face " << face;
    }
  }
}

TEST (t8_cmesh_cgns_reader, two_zones_hex)
{
  t8_cmesh_cgns_reader reader ("");
  t8_cmesh_t cmesh = t8_cgns_read_test_file ("test_cgns_two_zones.cgns", reader);
  ASSERT_NE (cmesh, nullptr);
  t8_cgns_check_two_zone_cmesh (cmesh, reader, false);
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_cgns_reader, two_structured_zones_hex)
{
  t8_cmesh_cgns_reader reader ("");
  t8_cmesh_t cmesh = t8_cgns_read_test_file ("test_cgns_two_structured_zones.cgns", reader);
  ASSERT_NE (cmesh, nullptr);
  t8_cgns_check_two_zone_cmesh (cmesh, reader, true);
  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_cgns_reader, mixed_2d)
{
  t8_cmesh_cgns_reader reader ("");
  t8_cmesh_t cmesh = t8_cgns_read_test_file ("test_cgns_mixed_2d.cgns", reader);
  ASSERT_NE (cmesh, nullptr);

  EXPECT_EQ (t8_cmesh_get_dimension (cmesh), 2);
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), 2);
  EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, 0), T8_ECLASS_QUAD);
  EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, 1), T8_ECLASS_TRIANGLE);
  /* The quad face 1 (x = 1, vertices (1,0) and (1,1)) is connected to the triangle face 0 (between its
   * vertices 1 = (1,1) and 2 = (1,0)). The triangle is the lower class, so its face is the main face and
   * its face corner 0 = (1,1) is the corner 1 of the quad face. */
  t8_cgns_expect_neighbor (cmesh, 0, 1, 1, 0, 1);
  t8_cgns_expect_neighbor (cmesh, 1, 0, 0, 1, 1);
  for (const int face : { 0, 2, 3 }) {
    EXPECT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 0, face)) << "quad face " << face;
  }
  for (const int face : { 1, 2 }) {
    EXPECT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 1, face)) << "triangle face " << face;
  }
  /* The z coordinate of all vertices is zero. */
  const double *vertices = t8_cmesh_get_tree_vertices (cmesh, 1);
  ASSERT_NE (vertices, nullptr);
  const double expected_vertices[9] = { 2, 0.5, 0, 1, 1, 0, 1, 0, 0 };
  for (int i = 0; i < 9; ++i) {
    EXPECT_DOUBLE_EQ (vertices[i], expected_vertices[i]) << "coordinate " << i;
  }
  /* Boundary conditions */
  ASSERT_EQ (reader.get_boundary_conditions ().size (), 1u);
  EXPECT_EQ (reader.get_boundary_conditions ()[0].name, "bottom");
  EXPECT_EQ (reader.get_boundary_conditions ()[0].num_faces, 1);
  EXPECT_EQ (t8_cmesh_get_boundary_condition (cmesh, 0, 2), "bottom");
  for (const int face : { 0, 1, 3 }) {
    EXPECT_EQ (t8_cmesh_get_boundary_condition (cmesh, 0, face), "") << "quad face " << face;
  }
  for (int face = 0; face < 3; ++face) {
    EXPECT_EQ (t8_cmesh_get_boundary_condition (cmesh, 1, face), "") << "triangle face " << face;
  }

  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_cgns_reader, trace_unstructured_periodic)
{
  t8_cmesh_cgns_reader reader ("");
  t8_cmesh_t cmesh = t8_cgns_read_test_file ("simpleElements_1-1.cgns", reader);
  ASSERT_NE (cmesh, nullptr);

  EXPECT_EQ (t8_cmesh_get_dimension (cmesh), 3);
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), 2);
  ASSERT_EQ (reader.get_zones ().size (), 1u);
  EXPECT_EQ (reader.get_zones ()[0].name, "dest0_block");
  EXPECT_EQ (reader.get_zones ()[0].num_trees, 2);
  EXPECT_EQ (reader.get_zones ()[0].num_vertices, 12);
  for (t8_locidx_t itree = 0; itree < 2; ++itree) {
    EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, itree), T8_ECLASS_HEX);
    EXPECT_EQ (t8_cmesh_get_cgns_zone (cmesh, itree), 0);
  }

  /* The two walls are the only domain boundaries. All other faces are connected to a neighbor,
   * most of them across the periodic interfaces of the one cell wide domain. */
  ASSERT_EQ (reader.get_boundary_conditions ().size (), 2u);
  for (size_t ibc = 0; ibc < 2; ++ibc) {
    const t8_cmesh_cgns_boundary_condition &bc = reader.get_boundary_conditions ()[ibc];
    EXPECT_EQ (bc.name, std::to_string (ibc + 1));
    EXPECT_EQ (bc.family, "Wall");
    EXPECT_EQ (bc.type, static_cast<int> (CGNS_ENUMV (FamilySpecified)));
    EXPECT_EQ (bc.zone, 0);
    EXPECT_EQ (bc.num_faces, 1);
  }
  int num_boundary_faces = 0;
  std::set<std::string> boundary_labels;
  for (t8_locidx_t itree = 0; itree < 2; ++itree) {
    for (int face = 0; face < 6; ++face) {
      const std::string label (t8_cmesh_get_boundary_condition (cmesh, itree, face));
      if (t8_cmesh_tree_face_is_boundary (cmesh, itree, face)) {
        ++num_boundary_faces;
        boundary_labels.insert (label);
        EXPECT_FALSE (label.empty ()) << "tree " << itree << " face " << face;
      }
      else {
        EXPECT_TRUE (label.empty ()) << "tree " << itree << " face " << face;
        int dual_face = -1, orientation = -1;
        const t8_locidx_t neighbor = t8_cmesh_get_face_neighbor (cmesh, itree, face, &dual_face, &orientation);
        EXPECT_GE (neighbor, 0) << "tree " << itree << " face " << face;
        EXPECT_LT (neighbor, 2) << "tree " << itree << " face " << face;
      }
    }
  }
  EXPECT_EQ (num_boundary_faces, 2);
  EXPECT_EQ (boundary_labels, (std::set<std::string> { "1", "2" }));

  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_cgns_reader, trace_structured_periodic)
{
  t8_cmesh_cgns_reader reader ("");
  t8_cmesh_t cmesh = t8_cgns_read_test_file ("box_9x9x9_peroidicY.cgns", reader);
  ASSERT_NE (cmesh, nullptr);

  /* The cells (i, j, k) of the 9 x 9 x 9 index grid are the trees i + 9 j + 81 k. */
  const int n = 9;
  const auto tree_id = [n] (int i, int j, int k) -> t8_locidx_t { return i + n * j + n * n * k; };
  EXPECT_EQ (t8_cmesh_get_dimension (cmesh), 3);
  ASSERT_EQ (t8_cmesh_get_num_trees (cmesh), n * n * n);
  ASSERT_EQ (reader.get_zones ().size (), 1u);
  EXPECT_EQ (reader.get_zones ()[0].name, "blk_0blk");
  EXPECT_EQ (reader.get_zones ()[0].num_trees, n * n * n);
  EXPECT_EQ (reader.get_zones ()[0].num_vertices, (n + 1) * (n + 1) * (n + 1));

  /* Boundary conditions: two walls and two "dead group" families, one on each boundary plane
   * in the i and k directions. The j direction is periodic. */
  ASSERT_EQ (reader.get_boundary_conditions ().size (), 4u);
  std::set<std::string> bc_names;
  for (size_t ibc = 0; ibc < 4; ++ibc) {
    const t8_cmesh_cgns_boundary_condition &bc = reader.get_boundary_conditions ()[ibc];
    bc_names.insert (bc.name);
    EXPECT_EQ (bc.family, ibc < 2 ? "Wall" : "dead group");
    EXPECT_EQ (bc.type, static_cast<int> (CGNS_ENUMV (FamilySpecified)));
    EXPECT_EQ (bc.zone, 0);
    EXPECT_EQ (bc.num_faces, n * n);
  }
  EXPECT_EQ (bc_names, (std::set<std::string> { "p1", "p2", "p3", "p4" }));

  for (int k = 0; k < n; ++k) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        const t8_locidx_t itree = tree_id (i, j, k);
        EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, itree), T8_ECLASS_HEX);
        EXPECT_EQ (t8_cmesh_get_cgns_zone (cmesh, itree), 0);
        /* Faces 0 and 1 are the i faces, 2 and 3 the j faces and 4 and 5 the k faces. */
        const bool expected_boundary[6] = { i == 0, i == n - 1, false, false, k == 0, k == n - 1 };
        for (int face = 0; face < 6; ++face) {
          EXPECT_EQ (t8_cmesh_tree_face_is_boundary (cmesh, itree, face), expected_boundary[face])
            << "tree " << itree << " face " << face;
          const std::string label (t8_cmesh_get_boundary_condition (cmesh, itree, face));
          EXPECT_EQ (!label.empty (), expected_boundary[face]) << "tree " << itree << " face " << face;
          if (expected_boundary[face]) {
            EXPECT_TRUE (bc_names.count (label) > 0) << "tree " << itree << " face " << face;
          }
        }
        /* Neighbors in the j direction, periodic at j = 0 and j = n - 1. */
        t8_cgns_expect_neighbor (cmesh, itree, 2, tree_id (i, (j + n - 1) % n, k), 3, 0);
        t8_cgns_expect_neighbor (cmesh, itree, 3, tree_id (i, (j + 1) % n, k), 2, 0);
        if (i > 0) {
          t8_cgns_expect_neighbor (cmesh, itree, 0, tree_id (i - 1, j, k), 1, 0);
        }
        if (k < n - 1) {
          t8_cgns_expect_neighbor (cmesh, itree, 5, tree_id (i, j, k + 1), 4, 0);
        }
      }
    }
  }

  t8_cmesh_destroy (&cmesh);
}

TEST (t8_cmesh_cgns_reader, missing_file)
{
  const std::string filename = (t8_test_data_dir / "test_cgns_does_not_exist.cgns").string ();
  t8_cmesh_t cmesh = t8_cmesh_from_cgns_file (filename.c_str (), sc_MPI_COMM_WORLD);
  EXPECT_EQ (cmesh, nullptr);
}
