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
#include <t8_vec.h>
#include <t8_element_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_examples.h>

#if T8_ENABLE_LESS_TESTS
#define MAXLEVEL 3
#else
#define MAXLEVEL 4
#endif

void
t8_element_centroid_by_vertex_coords (t8_forest_t forest, t8_eclass_scheme_c *ts, t8_locidx_t ltreeid,
                                      const t8_element_t *element, double *coordinates)
{
  double vertex_ref_coords[3], vertex_out_coords[3];
  t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);
  t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);

  /* Initialize the centroid with (0, 0, 0). */
  memset (coordinates, 0, 3 * sizeof (double));
  /* Get the number of corners of the element. */
  const int num_vertices = ts->t8_element_num_corners (element);
  for (int i_vertex = 0; i_vertex < num_vertices; i_vertex++) {
    /* For each corner, add its coordinates to the centroids coordinates. */

    /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
    ts->t8_element_vertex_reference_coords (element, i_vertex, vertex_ref_coords);
    /* Evaluate the geometry */
    t8_geometry_evaluate (cmesh, gtreeid, vertex_ref_coords, 1, vertex_out_coords);
    /* coordinates = coordinates + vertex_coords */
    t8_vec_axpy (vertex_out_coords, coordinates, 1);
  }
  /* Divide each coordinate by num_vertices */
  t8_vec_ax (coordinates, 1. / num_vertices);
}

void
t8_get_batch_coords_for_element_type (const t8_element_shape_t shape, double *batch_coords)
{
  const int num_vertices = t8_eclass_num_vertices[shape];
  const int elem_dim = t8_eclass_to_dimension[shape];
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    for (int dim = 0; dim < elem_dim; ++dim) {
      batch_coords[i_vertex * elem_dim + dim] = t8_element_corner_ref_coords[shape][i_vertex][dim];
    }
  }
}

template <typename T>
void
t8_write_message_by_dim (const char *message, const T *array, const int dim)
{
  std::ostringstream buffer;
  buffer.precision (6);
  buffer.setf (std::ios::fixed);
  buffer << message;
  for (int i_dim = 0; i_dim < dim; ++i_dim) {
    buffer << " " << array[i_dim];
  }
  t8_debugf ("%s\n", buffer.str ().c_str ());
}

/* To save an iterate, we will test the ref coords and centroid right in the adapt callback itself. 
 * Refines as long as the MAXLEVEL is not reached. */
int
t8_adapt_callback_with_test (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             t8_locidx_t lelement_id, t8_eclass_scheme_c *ts, const int is_family,
                             const int num_elements, t8_element_t *elements[])
{
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
  const t8_element_shape_t shape = ts->t8_element_shape (elements[0]);
  const int num_vertices = t8_eclass_num_vertices[shape];
  const int elem_dim = t8_eclass_to_dimension[shape];
  t8_get_batch_coords_for_element_type (shape, batch_coords);
  t8_debugf ("Testing for shape %s\n", t8_eclass_to_string[shape]);
  t8_debugf ("with num_vertices %i\n", num_vertices);
  t8_debugf ("and elem_dim %i\n", elem_dim);

  ts->t8_element_reference_coords (elements[0], batch_coords, num_vertices, tree_ref_coords_by_element_ref_coords);
  for (int i_vertex = 0; i_vertex < num_vertices; ++i_vertex) {
    ts->t8_element_vertex_reference_coords (elements[0], i_vertex, tree_ref_coords_by_vertex);
    t8_debugf ("Checking vertex %i\n", i_vertex);
    t8_write_message_by_dim ("with the batch coords:", batch_coords, elem_dim);
    t8_write_message_by_dim ("tree_ref_coords_by_vertex:", tree_ref_coords_by_vertex, elem_dim);
    t8_write_message_by_dim (
      "tree_ref_coords_by_element_ref_coords:", tree_ref_coords_by_element_ref_coords + i_vertex * elem_dim, elem_dim);
    for (int dim = 0; dim < elem_dim; ++dim) {
      EXPECT_NEAR (tree_ref_coords_by_vertex[dim], tree_ref_coords_by_element_ref_coords[i_vertex * elem_dim + dim],
                   2 * T8_PRECISION_EPS);
    }
  }
  /* compare results of the two different ways to compute an elements centroid */
  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid_by_element_ref_coords);
  t8_element_centroid_by_vertex_coords (forest_from, ts, which_tree, elements[0], centroid_by_vertices);
  t8_write_message_by_dim ("centroid_by_vertices:", centroid_by_vertices, 3);
  t8_write_message_by_dim ("centroid_by_element_ref_coords:", centroid_by_element_ref_coords, 3);
  for (int dim = 0; dim < T8_ECLASS_MAX_DIM; ++dim) {
    EXPECT_NEAR (centroid_by_vertices[dim], centroid_by_element_ref_coords[dim], 2 * T8_PRECISION_EPS);
  }

  /* refine if MAXLEVEL is not reached */
  if (ts->t8_element_level (elements[0]) >= MAXLEVEL) {
    return 0;
  }
  return 1;
}

class class_ref_coords: public testing::TestWithParam<t8_eclass> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
    forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_eclass_t eclass;
  t8_cmesh_t cmesh;
  t8_forest_t forest;
};

TEST_P (class_ref_coords, t8_check_elem_ref_coords)
{
  for (int level = 0; level < MAXLEVEL; ++level) {
    t8_forest_t forest_new;
    t8_forest_init (&forest_new);
    t8_forest_set_adapt (forest_new, forest, t8_adapt_callback_with_test, 0);
    t8_forest_set_partition (forest_new, forest, 0);
    t8_forest_commit (forest_new);
    forest = forest_new;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_ref_coords, class_ref_coords,
                          testing::Range (T8_ECLASS_VERTEX, T8_ECLASS_COUNT));