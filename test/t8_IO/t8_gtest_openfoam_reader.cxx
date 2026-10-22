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

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_openfoam/t8_openfoam_reader.hxx>
#include <test/t8_gtest_macros.hxx>
#include "t8_test_data_dir.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <vector>

/**
 * \file t8_gtest_openfoam_reader.cxx
 * In this file we test the OpenFOAM reader. The test cases are located in test/testfiles/openfoam.
 *
 * - single_<eclass>: A single cell cut out of the unit cube, see \ref unit_cube_cell_vertices.
 * - single_tet_uniform_owner: Like single_tet, but the owner list is written in the uniform encoding N{v}.
 * - error_*: Cases the reader has to reject (binary format, missing file, polyhedral cell, owner list with
 *   fewer values than announced).
 * - pitzDailySteadyCylinder_coarse: The pitzDaily tutorial with a cylinder cutout, meshed with half the
 *   resolution of the tutorial in each direction. The polyhedra created by snappyHexMesh are decomposed
 *   into pyramids and tets, so the mesh consists of hexes, prisms, pyramids and tets.
 *
 * The face connectivity computed by the reader is compared to the one computed by
 * \ref t8_cmesh_set_join_by_vertices from the tree vertices of the read cmesh. */

/** Tree vertices of a cmesh in the layout expected by \ref t8_cmesh_set_join_by_vertices. */
using t8_tree_vertices = std::array<t8_3D_vec, T8_ECLASS_MAX_CORNERS>;

/**
 * Read an OpenFOAM case from the test data directory.
 * \param [in] case_name  The name of the case directory in test/testfiles/openfoam.
 * \return                The read cmesh or NULL if the reader failed.
 */
static t8_cmesh_t
read_openfoam_test_case (const std::string& case_name)
{
  const std::filesystem::path foamfile = t8_test_data_dir / "openfoam" / case_name / "case.foam";
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_openfoam_reader reader (&cmesh, foamfile, sc_MPI_COMM_WORLD);
  reader.read ();
  return cmesh;
}

/**
 * The vertices of the single cell cases. All cells are cut out of the unit cube [0,1]^3, so all coordinates are 0 or 1.
 * \param [in] eclass  The eclass of the cell.
 * \return             The vertices of the cell in no particular order.
 */
static std::vector<t8_3D_vec>
unit_cube_cell_vertices (const t8_eclass_t eclass)
{
  /** The corners of the unit cube in t8code hex vertex order. */
  constexpr std::array<t8_3D_vec, 8> cube = { {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 1, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 1 },
    { 0, 1, 1 },
    { 1, 1, 1 },
  } };

  switch (eclass) {
  case T8_ECLASS_HEX:
    /* The unit cube itself. */
    return { cube[0], cube[1], cube[2], cube[3], cube[4], cube[5], cube[6], cube[7] };
  case T8_ECLASS_TET:
    /* The corner of the unit cube at the origin. */
    return { cube[0], cube[1], cube[2], cube[4] };
  case T8_ECLASS_PRISM:
    /* The half of the unit cube with x >= y. */
    return { cube[0], cube[1], cube[3], cube[4], cube[5], cube[7] };
  case T8_ECLASS_PYRAMID:
    /* The bottom face of the unit cube as base and the corner (0, 0, 1) as tip. */
    return { cube[0], cube[1], cube[2], cube[3], cube[4] };
  default:
    SC_ABORTF ("No single cell case for eclass %s.\n", t8_eclass_to_string[eclass]);
  }
}

/**
 * Check that a tree has the expected vertices.
 * Since the reader is free to choose the first vertex of a tree, the vertices are compared independent of their order.
 * \param [in] cmesh              A committed cmesh.
 * \param [in] itree              The local id of the tree to check.
 * \param [in] expected_vertices  The expected vertices of the tree in any order.
 */
static void
expect_tree_vertices (t8_cmesh_t cmesh, const t8_locidx_t itree, std::vector<t8_3D_vec> expected_vertices)
{
  const size_t num_vertices = t8_eclass_num_vertices[t8_cmesh_get_tree_class (cmesh, itree)];
  ASSERT_EQ (num_vertices, expected_vertices.size ());

  const double* vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
  std::vector<t8_3D_vec> tree_vertices (num_vertices);
  for (size_t ivertex = 0; ivertex < num_vertices; ++ivertex) {
    std::copy_n (vertices + 3 * ivertex, 3, tree_vertices[ivertex].begin ());
  }
  std::ranges::sort (tree_vertices);
  std::ranges::sort (expected_vertices);

  /* All coordinates are integers, which are read exactly, so we can compare without tolerance. */
  EXPECT_EQ (tree_vertices, expected_vertices);
}

/**
 * Check that the face connectivity of the cmesh matches the one computed by \ref t8_cmesh_set_join_by_vertices.
 * \param [in] cmesh  A committed, replicated cmesh.
 */
static void
expect_neighbors_match_join_by_vertices (t8_cmesh_t cmesh)
{
  const t8_locidx_t num_trees = t8_cmesh_get_num_local_trees (cmesh);

  /* Gather all eclasses and tree vertices. */
  std::vector<t8_eclass_t> eclasses (num_trees);
  std::vector<t8_tree_vertices> vertices (num_trees);
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    eclasses[itree] = t8_cmesh_get_tree_class (cmesh, itree);
    const double* tree_vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
    for (int ivertex = 0; ivertex < t8_eclass_num_vertices[eclasses[itree]]; ++ivertex) {
      std::copy_n (tree_vertices + 3 * ivertex, 3, vertices[itree][ivertex].begin ());
    }
  }

  /* Compute the face connectivity from the vertices. */
  int* connectivity = nullptr;
  t8_cmesh_set_join_by_vertices (nullptr, num_trees, eclasses.data (), vertices.data ()->data ()->data (),
                                 &connectivity, 1);

  /* Compare with the face connectivity of the cmesh. */
  for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
    for (int iface = 0; iface < t8_eclass_num_faces[eclasses[itree]]; ++iface) {
      const int* face_connectivity = connectivity + 3 * (T8_ECLASS_MAX_FACES * itree + iface);
      int dual_face = -1;
      int orientation = -1;
      const t8_locidx_t neighbor = t8_cmesh_get_face_neighbor (cmesh, itree, iface, &dual_face, &orientation);
      EXPECT_EQ (neighbor, face_connectivity[0]) << "Neighbor of tree " << itree << " face " << iface << " differs.";
      if (neighbor >= 0) {
        EXPECT_EQ (dual_face, face_connectivity[1])
          << "Dual face of tree " << itree << " face " << iface << " differs.";
        EXPECT_EQ (orientation, face_connectivity[2])
          << "Orientation of tree " << itree << " face " << iface << " differs.";
      }
    }
  }
  T8_FREE (connectivity);
}

/** Test fixture for the single cell cases, parameterized by the eclass of the cell. */
class t8_gtest_openfoam_reader_single_cell: public testing::TestWithParam<t8_eclass_t> {};

/**
 * Read the case single_<eclass>, which contains exactly one cell cut out of the unit cube.
 * This tests the reconstruction of each cell type on its own: The read cmesh must contain exactly one tree
 * with the expected eclass and the vertices of \ref unit_cube_cell_vertices (in any order).
 * Since the case has no internal faces, all faces of the tree must be boundary faces.
 */
TEST_P (t8_gtest_openfoam_reader_single_cell, unit_cube_cell)
{
  const t8_eclass_t eclass = GetParam ();
  std::string case_name = std::string ("single_") + t8_eclass_to_string[eclass];
  std::ranges::transform (case_name, case_name.begin (), ::tolower);

  t8_cmesh_t cmesh = read_openfoam_test_case (case_name);
  ASSERT_NE (cmesh, nullptr);

  ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), 1);
  EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, 0), eclass);
  expect_tree_vertices (cmesh, 0, unit_cube_cell_vertices (eclass));
  for (int iface = 0; iface < t8_eclass_num_faces[eclass]; ++iface) {
    EXPECT_TRUE (t8_cmesh_tree_face_is_boundary (cmesh, 0, iface));
  }

  t8_cmesh_destroy (&cmesh);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_openfoam_reader, t8_gtest_openfoam_reader_single_cell,
                          testing::Values (T8_ECLASS_HEX, T8_ECLASS_TET, T8_ECLASS_PRISM, T8_ECLASS_PYRAMID),
                          print_eclass);

/**
 * Read the case single_tet_uniform_owner, which is single_tet with the owner list written in the
 * uniform encoding "4{0}" instead of an expanded list. This tests the parsing of uniform label lists.
 * The read cmesh must contain exactly one tet with the vertices of \ref unit_cube_cell_vertices.
 */
TEST (t8_gtest_openfoam_reader, single_tet_uniform_owner)
{
  t8_cmesh_t cmesh = read_openfoam_test_case ("single_tet_uniform_owner");
  ASSERT_NE (cmesh, nullptr);

  ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), 1);
  EXPECT_EQ (t8_cmesh_get_tree_class (cmesh, 0), T8_ECLASS_TET);
  expect_tree_vertices (cmesh, 0, unit_cube_cell_vertices (T8_ECLASS_TET));

  t8_cmesh_destroy (&cmesh);
}

/** Test fixture for the cases the reader has to reject, parameterized by the name of the case. */
class t8_gtest_openfoam_reader_error: public testing::TestWithParam<std::string> {};

/**
 * Read a case the reader has to reject and check that it returns a NULL cmesh on all ranks instead of aborting.
 * All cases are based on single_hex:
 * - error_binary_format:     The files are declared as binary, which is not supported.
 * - error_missing_neighbour: The neighbour file does not exist.
 * - error_polyhedron:        One face of the hex is split into two triangles, so the cell is a polyhedron.
 * - error_truncated_owner:   The owner list announces 6 values but only contains 5.
 */
TEST_P (t8_gtest_openfoam_reader_error, returns_null)
{
  t8_cmesh_t cmesh = read_openfoam_test_case (GetParam ());
  EXPECT_EQ (cmesh, nullptr);
  if (cmesh != nullptr) {
    t8_cmesh_destroy (&cmesh);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_openfoam_reader, t8_gtest_openfoam_reader_error,
                          testing::Values ("error_binary_format", "error_missing_neighbour", "error_polyhedron",
                                           "error_truncated_owner"),
                          [] (const testing::TestParamInfo<std::string>& info) { return info.param; });

/**
 * Read the case pitzDailySteadyCylinder_coarse, which was generated by OpenFOAM (blockMesh, snappyHexMesh and a
 * decomposition of the resulting polyhedra into pyramids and tets). It is a realistic mesh with all four cell
 * types and highest cell ids which only appear in the neighbour file.
 * The number of trees per eclass must match the cell types reported by OpenFOAM's checkMesh and the face
 * connectivity must match \ref t8_cmesh_set_join_by_vertices.
 */
TEST (t8_gtest_openfoam_reader, pitz_daily)
{
  t8_cmesh_t cmesh = read_openfoam_test_case ("pitzDailySteadyCylinder_coarse");
  ASSERT_NE (cmesh, nullptr);

  /* Number of cells and cell types as reported by OpenFOAM's checkMesh. */
  ASSERT_EQ (t8_cmesh_get_num_local_trees (cmesh), 2284);
  std::array<int, T8_ECLASS_COUNT> num_trees_per_eclass {};
  for (t8_locidx_t itree = 0; itree < t8_cmesh_get_num_local_trees (cmesh); ++itree) {
    ++num_trees_per_eclass[t8_cmesh_get_tree_class (cmesh, itree)];
  }
  EXPECT_EQ (num_trees_per_eclass[T8_ECLASS_HEX], 1008);
  EXPECT_EQ (num_trees_per_eclass[T8_ECLASS_PRISM], 40);
  EXPECT_EQ (num_trees_per_eclass[T8_ECLASS_PYRAMID], 372);
  EXPECT_EQ (num_trees_per_eclass[T8_ECLASS_TET], 864);

  expect_neighbors_match_join_by_vertices (cmesh);
  t8_cmesh_destroy (&cmesh);
}
