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
/** \file t8_gtest_geometrical.cxx
 * Unit tests for the geometry of subelements. A single root element is refined uniformly to level 1, one element
 * is refined further and the hanging nodes are removed. For one subelement of a resulting transition cell, the
 * vertex coordinates (t8_forest_element_coordinate) and the face centroids (t8_forest_element_face_centroid) are
 * compared against values that are known by construction.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_adapt_callbacks.hxx>

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_subelement.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>
#include <t8_types/t8_vec.hxx>

#include <vector>

/** Build a level 1 forest with the subelement scheme on a single root element of class \a eclass, refine it with
 * \a adapt_callback and remove the hanging nodes. Then check the vertices and face centroids of the leaf \a ileaf.
 * \param [in] eclass             The class of the root element.
 * \param [in] adapt_callback     The callback refining the elements that cause the hanging nodes.
 * \param [in] ileaf              The index of the leaf to check in the (only) tree after adapt. 
 *                                Should be a subelement index.
 * \param [in] expected_vertex    The expected coordinates of all vertices of the leaf.
 * \param [in] expected_centroid  The expected coordinates of all face centroids of the leaf.
 */
static void
check_subelement_geometry (const t8_eclass_t eclass, t8_forest_adapt_t adapt_callback, const t8_locidx_t ileaf,
                           const std::vector<t8_3D_vec> &expected_vertex,
                           const std::vector<t8_3D_vec> &expected_centroid)
{
  int mpisize;
  int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  if (!(mpisize == 1)) {
    GTEST_SKIP () << "Skipping test if run in parallel.";
  }

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_from_class (cmesh, eclass, sc_MPI_COMM_WORLD);
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_subelement (), 1, 0, sc_MPI_COMM_WORLD);
  forest = t8_forest_new_adapt (forest, adapt_callback, 0, 0, NULL);
  forest = t8_forest_remove_hanging_nodes (forest);
  EXPECT_TRUE (t8_forest_has_subelements (forest));

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, 0);
  EXPECT_LT (ileaf, t8_forest_get_tree_num_leaf_elements (forest, 0));
  const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, 0, ileaf);
  EXPECT_TRUE (t8_element_is_subelement (scheme, tree_class, elem));

  /* Vertices. */
  const int num_corners = scheme->element_get_num_corners (tree_class, elem);
  EXPECT_EQ (num_corners, static_cast<int> (expected_vertex.size ()));
  for (int icorner = 0; icorner < num_corners; ++icorner) {
    t8_3D_vec vertex {};
    t8_forest_element_coordinate (forest, 0, elem, icorner, vertex.data ());
    for (int d = 0; d < 3; ++d) {
      EXPECT_NEAR (vertex[d], expected_vertex[icorner][d], T8_PRECISION_SQRT_EPS)
        << "vertex " << icorner << ", dim " << d;
    }
  }

  /* Face centroids. */
  const int num_faces = scheme->element_get_num_faces (tree_class, elem);
  EXPECT_EQ (num_faces, static_cast<int> (expected_centroid.size ()));
  for (int iface = 0; iface < num_faces; ++iface) {
    t8_3D_vec centroid {};
    t8_forest_element_face_centroid (forest, 0, elem, iface, centroid.data ());
    for (int d = 0; d < 3; ++d) {
      EXPECT_NEAR (centroid[d], expected_centroid[iface][d], T8_PRECISION_SQRT_EPS)
        << "face " << iface << ", dim " << d;
    }
  }

  t8_forest_unref (&forest);
}

/** Unit square, level 1. Refining element 2 (upper left) makes only the upper face of the lower left quad hanging
 * (type 1, 5 subelements). Its subelement 1 (leaf 1) is the left half of the hanging face with vertices:
 * centre (0.25, 0.25), upper left corner (0, 0.5) and face midpoint (0.25, 0.5). */
TEST (t8_gtest_subelement_geometry, quad_vertices_and_face_centroids)
{
  check_subelement_geometry (T8_ECLASS_QUAD, refine_every_nth_element_callback<4, 2>, 1,
                             { { 0.25, 0.25, 0.0 }, { 0.0, 0.5, 0.0 }, { 0.25, 0.5, 0.0 } },
                             { { 0.125, 0.5, 0.0 }, { 0.25, 0.375, 0.0 }, { 0.125, 0.375, 0.0 } });
}

/** Reference triangle (0,0), (1,0), (1,1), level 1. Refining element 2 (middle child) makes only f0 of the corner
 * child (0,0), (0.5,0), (0.5,0.5) hanging (type 4, 2 subelements). Its subelement 0 (leaf 0) has the vertices:
 * midpoint of f0 (0.5, 0.25), v1 (0.5, 0) and v0 (0, 0). */
TEST (t8_gtest_subelement_geometry, tri_vertices_and_face_centroids)
{
  check_subelement_geometry (T8_ECLASS_TRIANGLE, refine_every_nth_element_callback<4, 2>, 0,
                             { { 0.5, 0.25, 0.0 }, { 0.5, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
                             { { 0.25, 0.0, 0.0 }, { 0.25, 0.125, 0.0 }, { 0.5, 0.125, 0.0 } });
}
