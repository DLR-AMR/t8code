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
#include <t8_eclass.h>
#include <t8_types/t8_vec.h>
#include <t8_element.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <test/t8_gtest_macros.hxx>

#if T8CODE_TEST_LEVEL >= 1
#define MAX_LEVEL_REF_COORD_TEST 3
#else
#define MAX_LEVEL_REF_COORD_TEST 4
#endif

/**
 * Writes the contents of an array into a string.
 * @tparam T The type of the array.
 * \param [in,out] message The message which is followed by the contents of the array.
 * \param [in,out] array The array to write.
 * \param [in,out] dim The dimension of the array.
 * \return The string containing the message and the array.
 */
template <typename T>
std::string
t8_write_message_by_dim (const char *message, const T *array, const int dim)
{
  std::ostringstream buffer;
  buffer.precision (10);
  buffer.setf (std::ios::fixed);
  buffer << message;
  for (int i_dim = 0; i_dim < dim; ++i_dim) {
    buffer << " " << array[i_dim];
  }
  return buffer.str ();
}

/**
 * Computes the centroid of an element by computing the coordinates of the vertices and computing the mean of them.
 * \param [in] forest The forest.
 * \param [in] ltreeid The local tree id.
 * \param [in] element The element.
 * \param [out] coordinates The coordinates of the centroid.
 */
void
t8_element_centroid_by_vertex_coords (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                      double *coordinates)
{
  double vertex_ref_coords[3], vertex_out_coords[3];
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);

  /* Initialize the centroid with (0, 0, 0). */
  std::fill_n (coordinates, 3, 0);
  /* Get the number of corners of the element. */
  const int num_vertices = scheme->element_get_num_corners (tree_class, element);
  for (int i_vertex = 0; i_vertex < num_vertices; i_vertex++) {
    /* For each corner, add its coordinates to the centroids coordinates. */

    /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
    scheme->element_get_vertex_reference_coords (tree_class, element, i_vertex, vertex_ref_coords);
    /* Evaluate the geometry */
    t8_geometry_evaluate (cmesh, gtreeid, vertex_ref_coords, 1, vertex_out_coords);
    /* coordinates = coordinates + vertex_coords */
    t8_axpy_c_interface (vertex_out_coords, coordinates, 1);
  }
  /* Divide each coordinate by num_vertices */
  t8_ax_c_interface (coordinates, 1. / num_vertices);
}

/**
 * Computes the reference coordinates of the vertices of an element.
 * \param [in] shape The element shape.
 * \param [out] batch_coords The batch coordinates of the vertices for the element shape.
 */
void
t8_get_batch_coords_for_element_type (const t8_element_shape_t shape, double *batch_coords)
{
  const int num_vertices = t8_eclass_num_vertices[shape];
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
      batch_coords[i_vertex * T8_ECLASS_MAX_DIM + dim] = t8_element_corner_ref_coords[shape][i_vertex][dim];
    }
  }
}

/**
 * Generates additional info for the reference coordinates test.
 * \param [in] shape The shape of the element.
 * \param [in] i_vertex The vertex index.
 * \param [in] elem_dim The dimension of the element.
 * \param [in] batch_coords The input batch coordinates.
 * \param [in] tree_ref_coords_by_vertex The reference coordinates of the vertex computed by \ref t8_element_vertex_reference_coords.
 * \param [in] tree_ref_coords_by_element_ref_coords The reference coordinates of the vertex computed by \ref t8_element_reference_coords.
 * \return The additional info.
 */
std::string
t8_generate_additional_info_ref_coords (const t8_element_shape_t shape, const int i_vertex, const int elem_dim,
                                        const double *batch_coords, const double *tree_ref_coords_by_vertex,
                                        const double *tree_ref_coords_by_element_ref_coords)
{
  std::ostringstream add_info;
  add_info << "Test failed for element shape " << t8_eclass_to_string[shape];
  add_info << " on vertex " << i_vertex << std::endl;
  add_info << t8_write_message_by_dim ("with the batch coords:", batch_coords, elem_dim) << std::endl;
  add_info << t8_write_message_by_dim ("tree_ref_coords_by_vertex:", tree_ref_coords_by_vertex, elem_dim) << std::endl;
  add_info << t8_write_message_by_dim (
    "tree_ref_coords_by_element_ref_coords:", tree_ref_coords_by_element_ref_coords + i_vertex * elem_dim, elem_dim);
  return add_info.str ();
}

/**
 * Generates additional info for the centroid test.
 * \param [in] shape The shape of the element.
 * \param [in] centroid_by_vertices The centroid computed by \ref t8_element_vertex_reference_coords -> \ref t8_geometry_evaluate -> mean of all results.
 * \param [in] centroid_by_element_ref_coords The centroid computed by \ref t8_forest_element_centroid (uses \ref t8_element_reference_coords) -> \ref t8_geometry_evaluate.
 * \return The additional info.
 */
std::string
t8_generate_additional_info_centroid (const t8_element_shape_t shape, const double *centroid_by_vertices,
                                      const double *centroid_by_element_ref_coords)
{
  std::ostringstream add_info;
  add_info << "Test failed for element shape " << t8_eclass_to_string[shape] << std::endl;
  add_info << t8_write_message_by_dim ("centroid_by_vertices:", centroid_by_vertices, 3) << std::endl;
  add_info << t8_write_message_by_dim ("centroid_by_element_ref_coords:", centroid_by_element_ref_coords, 3);
  return add_info.str ();
}

/**
 * Compares two arrays of doubles with a given tolerance.
 * \param [in] array1 The first array to compare.
 * \param [in] array2 The second array to compare.
 * \param [in] dim    The dimension of the arrays.
 * \param [in] tol    The tolerance.
 * \return True if the arrays are equal, false otherwise.
 */
const bool
t8_compare_arrays (const double *array1, const double *array2, const int dim, const double tol)
{
  for (int i_dim = 0; i_dim < dim; ++i_dim) {
    if (std::abs (array1[i_dim] - array2[i_dim]) > tol) {
      return false;
    }
  }
  return true;
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

  double tree_ref_coords_by_vertex
    [3]; /** reference coordinates of the element vertices computed by \ref t8_element_vertex_reference_coords */
  double tree_ref_coords_by_element_ref_coords
    [T8_ECLASS_MAX_CORNERS
     * T8_ECLASS_MAX_DIM]; /** reference coordinates of the element vertices computed by \ref t8_element_reference_coords */
  double centroid_by_vertices
    [3]; /** centroid computed via \ref t8_element_vertex_reference_coords -> \ref t8_geometry_evaluate -> mean of all results */
  double centroid_by_element_ref_coords
    [3]; /** centroid computed via \ref t8_forest_element_centroid (uses \ref t8_element_reference_coords) -> \ref t8_geometry_evaluate */
  double batch_coords
    [T8_ECLASS_MAX_CORNERS
     * T8_ECLASS_MAX_DIM]; /** reference coordinates of the element vertices computed by \ref t8_get_batch_coords_for_element_type */

  /* compare results of the two different way to obtain tree ref coords */
  const t8_element_shape_t shape = scheme->element_get_shape (tree_class, element);
  const int num_vertices = t8_eclass_num_vertices[shape];
  const int elem_dim = t8_eclass_to_dimension[shape];
  t8_get_batch_coords_for_element_type (shape, batch_coords);

  scheme->element_get_reference_coords (tree_class, element, batch_coords, num_vertices,
                                        tree_ref_coords_by_element_ref_coords);
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    scheme->element_get_vertex_reference_coords (tree_class, element, i_vertex, tree_ref_coords_by_vertex);
    EXPECT_TRUE (t8_compare_arrays (tree_ref_coords_by_vertex,
                                    tree_ref_coords_by_element_ref_coords + i_vertex * elem_dim, elem_dim,
                                    2 * T8_PRECISION_EPS))
      << t8_generate_additional_info_ref_coords (shape, i_vertex, elem_dim, batch_coords, tree_ref_coords_by_vertex,
                                                 tree_ref_coords_by_element_ref_coords);
  }
  /* Compare results of the two different ways to compute an elements centroid */
  t8_forest_element_centroid (forest, ltree_id, element, centroid_by_element_ref_coords);
  t8_element_centroid_by_vertex_coords (forest, ltree_id, element, centroid_by_vertices);
  EXPECT_TRUE (
    t8_compare_arrays (centroid_by_vertices, centroid_by_element_ref_coords, T8_ECLASS_MAX_DIM, 2 * T8_PRECISION_EPS))
    << t8_generate_additional_info_centroid (shape, centroid_by_vertices, centroid_by_element_ref_coords);
}

class class_ref_coords: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
 protected:
  void
  SetUp () override
  {
    const std::tuple<t8_eclass, int> params = GetParam ();
    const t8_eclass_t eclass = std::get<0> (params);
    const int level = std::get<1> (params);
    t8_cmesh_t cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
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
    for (ielement = 0; ielement < t8_forest_get_tree_num_elements (forest, itree); ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      t8_test_coords (forest, itree, element);
    }
  }
  /* Increase cmesh ref counter to not loose it during t8_forest_unref */
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_ref_coords, class_ref_coords,
                          testing::Combine (AllEclasses, testing::Range (0, MAX_LEVEL_REF_COORD_TEST + 1)));
