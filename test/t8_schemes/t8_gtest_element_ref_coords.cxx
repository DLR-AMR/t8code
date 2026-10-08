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

/** \file t8_gtest_element_ref_coords.cxx
* Provide tests to check the functionality of the computation of
* element reference coordinates and element centroids.
*/

#include <gtest/gtest.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_types/t8_vec.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_scheme.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_geometrical.hxx>
#include <t8_geometry/t8_geometry.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <sstream>
#include <vector>

#if T8_TEST_LEVEL_INT >= 1
#define MAX_LEVEL_REF_COORD_TEST 3
#else
#define MAX_LEVEL_REF_COORD_TEST 4
#endif

/**
 * Writes the coordinates of a point into a string.
 * \param [in] message The message which is followed by the coordinates of the point.
 * \param [in] point   The point to write.
 * \return The string containing the message and the coordinates.
 */
std::string
t8_write_message_by_point (const char *message, const t8_3D_vec &point)
{
  std::ostringstream buffer;
  buffer.precision (10);
  buffer.setf (std::ios::fixed);
  buffer << message;
  for (const double coord : point) {
    buffer << " " << coord;
  }
  return buffer.str ();
}

/**
 * Computes the centroid of an element by computing the coordinates of the vertices and computing the mean of them.
 * \param [in] forest The forest.
 * \param [in] ltreeid The local tree id.
 * \param [in] element The element.
 * \return The coordinates of the centroid.
 */
t8_3D_vec
t8_element_centroid_by_vertex_coords (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element)
{
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);

  t8_3D_vec centroid = { 0, 0, 0 };
  /* Get the number of corners of the element. */
  const int num_vertices = scheme->element_get_num_corners (tree_class, element);
  for (int i_vertex = 0; i_vertex < num_vertices; i_vertex++) {
    /* For each corner, add its coordinates to the centroids coordinates. */

    /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
    const t8_3D_vec vertex_ref_coords = scheme->element_get_vertex_reference_coords (tree_class, element, i_vertex);
    /* Evaluate the geometry */
    const t8_3D_vec vertex_out_coords = t8_geometry_evaluate (cmesh, gtreeid, vertex_ref_coords);
    /* centroid = centroid + vertex_coords */
    t8_axpy (vertex_out_coords, centroid, 1);
  }
  /* Divide each coordinate by num_vertices */
  t8_ax (centroid, 1. / num_vertices);
  return centroid;
}

/**
 * Computes the reference coordinates of the vertices of an element.
 * \param [in] shape The element shape.
 * \return The reference coordinates of the vertices of the element shape.
 */
std::vector<t8_3D_vec>
t8_get_batch_coords_for_element_type (const t8_element_shape_t shape)
{
  const t8_3D_vec *corner_ref_coords = t8_element_corner_ref_coords[shape];
  return std::vector<t8_3D_vec> (corner_ref_coords, corner_ref_coords + t8_eclass_num_vertices[shape]);
}

/**
 * Generates additional info for the reference coordinates test.
 * \param [in] shape The shape of the element.
 * \param [in] i_vertex The vertex index.
 * \param [in] batch_coords The input batch coordinates of the vertex.
 * \param [in] tree_ref_coords_by_vertex The reference coordinates of the vertex computed by \ref t8_scheme::element_get_vertex_reference_coords.
 * \param [in] tree_ref_coords_by_element_ref_coords The reference coordinates of the vertex computed by \ref t8_scheme::element_get_reference_coords.
 * \return The additional info.
 */
std::string
t8_generate_additional_info_ref_coords (const t8_element_shape_t shape, const int i_vertex,
                                        const t8_3D_vec &batch_coords, const t8_3D_vec &tree_ref_coords_by_vertex,
                                        const t8_3D_vec &tree_ref_coords_by_element_ref_coords)
{
  std::ostringstream add_info;
  add_info << "Test failed for element shape " << t8_eclass_to_string[shape];
  add_info << " on vertex " << i_vertex << std::endl;
  add_info << t8_write_message_by_point ("with the batch coords:", batch_coords) << std::endl;
  add_info << t8_write_message_by_point ("tree_ref_coords_by_vertex:", tree_ref_coords_by_vertex) << std::endl;
  add_info << t8_write_message_by_point ("tree_ref_coords_by_element_ref_coords:",
                                         tree_ref_coords_by_element_ref_coords);
  return add_info.str ();
}

/**
 * Generates additional info for the centroid test.
 * \param [in] shape The shape of the element.
 * \param [in] centroid_by_vertices The centroid computed by \ref t8_scheme::element_get_vertex_reference_coords -> \ref t8_geometry_evaluate -> mean of all results.
 * \param [in] centroid_by_element_ref_coords The centroid computed by \ref t8_forest_element_centroid (uses \ref t8_scheme::element_get_reference_coords) -> \ref t8_geometry_evaluate.
 * \return The additional info.
 */
std::string
t8_generate_additional_info_centroid (const t8_element_shape_t shape, const t8_3D_vec &centroid_by_vertices,
                                      const t8_3D_vec &centroid_by_element_ref_coords)
{
  std::ostringstream add_info;
  add_info << "Test failed for element shape " << t8_eclass_to_string[shape] << std::endl;
  add_info << t8_write_message_by_point ("centroid_by_vertices:", centroid_by_vertices) << std::endl;
  add_info << t8_write_message_by_point ("centroid_by_element_ref_coords:", centroid_by_element_ref_coords);
  return add_info.str ();
}

/**
 * Tests the reference coordinates of an element.
 * \param [in] forest The forest.
 * \param [in] ltree_id The local tree id.
 * \param [in] element The element.
 */
void
t8_test_coords (const t8_forest_t forest, const t8_locidx_t ltree_id, const t8_element_t *element)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltree_id);
  const t8_element_shape_t shape = scheme->element_get_shape (tree_class, element);
  const int num_vertices = t8_eclass_num_vertices[shape];
  const int elem_dim = t8_eclass_to_dimension[shape];

  /** Reference coordinates of the element vertices. Components exceeding the dimension of the element have to be
   * ignored, hence we set them to an arbitrary value. */
  std::vector<t8_3D_vec> batch_coords = t8_get_batch_coords_for_element_type (shape);
  for (t8_3D_vec &coords : batch_coords) {
    std::fill (coords.begin () + elem_dim, coords.end (), 0.123);
  }

  /** Reference coordinates of the element vertices computed by \ref t8_scheme::element_get_reference_coords.
   * Initialized with -1 to check that the components exceeding the dimension are set to NaN in debug mode. */
  std::vector<t8_3D_vec> tree_ref_coords_by_element_ref_coords (num_vertices, t8_3D_vec { -1, -1, -1 });
  scheme->element_get_reference_coords (tree_class, element, batch_coords, tree_ref_coords_by_element_ref_coords);

  /** The same reference coordinates computed via the C interface, which stores each point as 3 doubles. */
  std::vector<double> c_tree_ref_coords (3 * num_vertices, -1);
  {
    std::vector<double> c_batch_coords (3 * num_vertices);
    t8_3D_vecs_to_doubles (batch_coords, c_batch_coords.data ());
    t8_element_get_reference_coords (scheme, tree_class, element, c_batch_coords.data (), num_vertices,
                                     c_tree_ref_coords.data ());
  }
  const std::vector<t8_3D_vec> tree_ref_coords_by_c_interface
    = t8_3D_vecs_from_doubles (c_tree_ref_coords.data (), num_vertices);

  /* compare results of the different ways to obtain tree ref coords */
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    /** Reference coordinates of the vertex computed by \ref t8_scheme::element_get_vertex_reference_coords. */
    const t8_3D_vec tree_ref_coords_by_vertex
      = scheme->element_get_vertex_reference_coords (tree_class, element, i_vertex);
    const t8_3D_vec &tree_ref_coords_by_batch = tree_ref_coords_by_element_ref_coords[i_vertex];
    /* Only the first elem_dim components are part of the result of a batch. */
    EXPECT_TRUE (t8_eq (tree_ref_coords_by_vertex | std::views::take (elem_dim),
                        tree_ref_coords_by_batch | std::views::take (elem_dim), 2 * T8_PRECISION_EPS))
      << t8_generate_additional_info_ref_coords (shape, i_vertex, batch_coords[i_vertex], tree_ref_coords_by_vertex,
                                                 tree_ref_coords_by_batch);
    EXPECT_TRUE (std::ranges::equal (tree_ref_coords_by_c_interface[i_vertex] | std::views::take (elem_dim),
                                     tree_ref_coords_by_batch | std::views::take (elem_dim)));
#if T8_ENABLE_DEBUG
    /* In debug mode, components exceeding the dimension of the element are NaN. */
    for (int i_dim = elem_dim; i_dim < T8_ECLASS_MAX_DIM; ++i_dim) {
      EXPECT_TRUE (std::isnan (tree_ref_coords_by_vertex[i_dim]));
      EXPECT_TRUE (std::isnan (tree_ref_coords_by_c_interface[i_vertex][i_dim]));
      EXPECT_TRUE (std::isnan (tree_ref_coords_by_batch[i_dim]));
    }
#endif
  }
  /* Compare results of the two different ways to compute an elements centroid */
  t8_3D_vec centroid_by_element_ref_coords;
  t8_forest_element_centroid (forest, ltree_id, element, centroid_by_element_ref_coords.data ());
  const t8_3D_vec centroid_by_vertices = t8_element_centroid_by_vertex_coords (forest, ltree_id, element);
  EXPECT_TRUE (t8_eq (centroid_by_vertices, centroid_by_element_ref_coords, 2 * T8_PRECISION_EPS))
    << t8_generate_additional_info_centroid (shape, centroid_by_vertices, centroid_by_element_ref_coords);
}

struct class_ref_coords: public testing::TestWithParam<std::tuple<t8_eclass_t, int>>
{
 protected:
  void
  SetUp () override
  {
    const std::tuple<t8_eclass, int> params = GetParam ();
    const t8_eclass_t eclass = std::get<0> (params);
    const int level = std::get<1> (params);
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_new_from_class (cmesh, eclass, sc_MPI_COMM_WORLD);
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), level, 0, sc_MPI_COMM_WORLD);
    t8_forest_init (&forest_partition);
    t8_forest_set_partition (forest_partition, forest, 0);
    t8_forest_commit (forest_partition);
    forest = forest_partition;
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_forest_t forest, forest_partition;
};

TEST_P (class_ref_coords, t8_check_elem_ref_coords)
{
  t8_locidx_t itree, ielement;
  /* Check the reference coordinates of each element in each tree */
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    for (ielement = 0; ielement < t8_forest_get_tree_num_leaf_elements (forest, itree); ielement++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      t8_test_coords (forest, itree, element);
    }
  }
  /* Increase cmesh ref counter to not loose it during t8_forest_unref */
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_ref_coords, class_ref_coords,
                          testing::Combine (AllEclasses, testing::Range (0, MAX_LEVEL_REF_COORD_TEST + 1)));
