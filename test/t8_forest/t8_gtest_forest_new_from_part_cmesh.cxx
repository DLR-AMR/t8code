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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_private.h>
#include "test/t8_cmesh_generator/t8_cmesh_example_sets.hxx"
#include <test/t8_gtest_macros.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
TEST (t8_forest_new_from_partition, all_elems_on_root)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  /* create a line cmesh */
  int mpirank, mpisize;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  sc_MPI_Comm_size (comm, &mpisize);
  sc_MPI_Comm_rank (comm, &mpirank);
  int face_knowledge = 3;
  int num_trees = 2;
  if (mpirank == 0) {
    t8_gloidx_t first = num_trees * mpirank;
    t8_gloidx_t second = num_trees * mpirank + 1;
    t8_gloidx_t third = num_trees * mpirank + 2;
    double vertices[9] = { (double) first, 0, 0, (double) second, 0, 0, (double) third, 0, 0 };

    t8_cmesh_set_tree_class (cmesh, first, T8_ECLASS_LINE);
    t8_cmesh_set_tree_class (cmesh, second, T8_ECLASS_LINE);
    t8_cmesh_set_tree_vertices (cmesh, first, vertices, 2);
    t8_cmesh_set_tree_vertices (cmesh, second, vertices + 3, 2);
    //    t8_cmesh_set_join (cmesh, first, second, 1, 0, 0);
    t8_cmesh_set_partition_range (cmesh, face_knowledge, 0, 1);
  }
  else {
    t8_cmesh_set_dimension (cmesh, 1);
    t8_cmesh_set_partition_range (cmesh, face_knowledge, 2, 1);
  }
  t8_cmesh_commit (cmesh, comm);
  EXPECT_EQ (t8_cmesh_get_num_local_trees (cmesh), (mpirank == 0) ? 2 : 0);
  EXPECT_EQ (t8_cmesh_get_num_trees (cmesh), 2);
  t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
  int level = 3;
  t8_forest_t forest = t8_forest_new_uniform_with_cmesh_partition (cmesh, scheme, level, 0, comm);
  EXPECT_EQ (t8_forest_get_local_num_elements (forest), (mpirank == 0) ? 2 << level : 0);
  EXPECT_EQ (t8_forest_get_global_num_elements (forest), 2 << level);
  t8_forest_unref (&forest);
}

TEST (t8_forest_new_from_partition, alternating_two_one)
{
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  /* create a line cmesh */
  int mpirank, mpisize;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  sc_MPI_Comm_size (comm, &mpisize);
  sc_MPI_Comm_rank (comm, &mpirank);
  int face_knowledge = 3;
  int num_trees_before = 3 * (mpirank / 2) + 2 * (mpirank % 2);
  t8_debugf ("num_trees_before:%i\n", num_trees_before);
  t8_gloidx_t first = num_trees_before;
  t8_gloidx_t second = num_trees_before + 1;
  t8_gloidx_t third = num_trees_before + 2;
  double vertices[9] = { (double) first, 0, 0, (double) second, 0, 0, (double) third, 0, 0 };
  t8_cmesh_set_tree_class (cmesh, first, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, first, vertices, 2);
  if (mpirank % 2 == 0) {
    t8_cmesh_set_tree_class (cmesh, second, T8_ECLASS_LINE);
    t8_cmesh_set_tree_vertices (cmesh, second, vertices + 3, 2);
  }
  t8_cmesh_set_partition_range (cmesh, face_knowledge, num_trees_before, num_trees_before + 1 - mpirank % 2);
  t8_cmesh_commit (cmesh, comm);
  int expected_local_trees = 2 - mpirank % 2;
  int expected_global_trees = 3 * (mpisize / 2) + 2 * (mpisize % 2);
  EXPECT_EQ (t8_cmesh_get_num_local_trees (cmesh), expected_local_trees);
  EXPECT_EQ (t8_cmesh_get_num_trees (cmesh), expected_global_trees);
  t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
  int level = 3;
  t8_forest_t forest = t8_forest_new_uniform_with_cmesh_partition (cmesh, scheme, level, 0, comm);
  EXPECT_EQ (t8_forest_get_local_num_elements (forest), expected_local_trees << level);
  EXPECT_EQ (t8_forest_get_global_num_elements (forest), expected_global_trees << level);

  t8_eclass_scheme_c *line_scheme = t8_forest_get_eclass_scheme (forest, T8_ECLASS_LINE);

  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    t8_tree_t tree = t8_forest_get_tree (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < t8_forest_get_tree_element_count (tree); ielement++) {
      const t8_element_t *element = t8_forest_get_tree_element (tree, ielement);
      EXPECT_EQ (line_scheme->t8_element_level (element), level);
    }
  }
  t8_forest_unref (&forest);
}

TEST (t8_forest_new_from_partition, empty_test)
{
  /** TODO:
   * hybrid
   * parameterized above over eclass
   * Look in cmesh examples, if there are any than can be constructed partitioned imbalanced.
  */
}

TEST (t8_forest_new_from_partition, one_two_quad)
{
  /* Current process */
  int rank, size;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);

  /* Geometry */
  int degree = 2;
  std::vector<std::vector<double>> cells;
  cells.push_back (std::vector<double> { -1, 1.2, 0,   0, 1.2,  0,   -1, 3,    0, 0, 3,    0,   -1, 2.1,
                                         0,  0,   2.1, 0, -0.5, 1.2, 0,  -0.5, 3, 0, -0.5, 2.1, 0 });
  cells.push_back (std::vector<double> { 0, 1.2, 0,   2, 1.5, 0,    0, 3, 0, 2, 3,   0,   0, 2.1,
                                         0, 1.7, 2.1, 0, 1,   1.55, 0, 1, 3, 0, 0.9, 2.1, 0 });
  cells.push_back (
    std::vector<double> { 2, 1.5, 0, 4, 0, 0, 2, 3, 0, 4, 3, 0, 1.7, 2.1, 0, 4, 1.5, 0, 3, 1, 0, 3, 3, 0, 3, 2, 0 });

  /* Create a partitioned cmesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  if (size != 2)
    SC_ABORT ("Program must be run on two MPI processes.\n");
  else {
    if (rank == 0) {
      t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_vertices (cmesh, 0, cells[0].data (), cells[0].size ());
      t8_cmesh_set_tree_vertices (cmesh, 1, cells[1].data (), cells[1].size ());
      t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
      t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
      t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
      t8_cmesh_set_partition_range (cmesh, 3, 0, 0);
    }
    else if (rank == 1) {
      for (size_t itree = 0; itree < 3; ++itree) {
        t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
        t8_cmesh_set_tree_vertices (cmesh, itree, cells[itree].data (), cells[itree].size ());
        t8_cmesh_set_attribute (cmesh, itree, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree,
                                sizeof (int), 1);
      }
      t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
      t8_cmesh_set_partition_range (cmesh, 3, 1, 2);
    }
  }
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh, 2);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest
    = t8_forest_new_uniform_with_cmesh_partition (cmesh, t8_scheme_new_default_cxx (), 0, 0, sc_MPI_COMM_WORLD);
  // The following line still works, because the partition "by chance" matches the partition, that t8code expects from
  // an equidistributed uniform mesh
  //   t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 1, sc_MPI_COMM_WORLD);
  t8_forest_unref (&forest);
}

TEST (t8_forest_new_from_partition, two_one_quad)
{
  /* Current process */
  int rank, size;
  sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
  sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &rank);

  /* Geometry */
  int degree = 2;
  std::vector<std::vector<double>> cells;
  cells.push_back (std::vector<double> { -1, 1.2, 0,   0, 1.2,  0,   -1, 3,    0, 0, 3,    0,   -1, 2.1,
                                         0,  0,   2.1, 0, -0.5, 1.2, 0,  -0.5, 3, 0, -0.5, 2.1, 0 });
  cells.push_back (std::vector<double> { 0, 1.2, 0,   2, 1.5, 0,    0, 3, 0, 2, 3,   0,   0, 2.1,
                                         0, 1.7, 2.1, 0, 1,   1.55, 0, 1, 3, 0, 0.9, 2.1, 0 });
  cells.push_back (
    std::vector<double> { 2, 1.5, 0, 4, 0, 0, 2, 3, 0, 4, 3, 0, 1.7, 2.1, 0, 4, 1.5, 0, 3, 1, 0, 3, 3, 0, 3, 2, 0 });

  /* Create a partitioned cmesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  if (size != 2)
    SC_ABORT ("Program must be run on two MPI processes.\n");
  else {
    if (rank == 0) {
      for (size_t itree = 0; itree < 3; ++itree) {
        t8_cmesh_set_tree_class (cmesh, itree, T8_ECLASS_QUAD);
        t8_cmesh_set_tree_vertices (cmesh, itree, cells[itree].data (), cells[itree].size ());
        t8_cmesh_set_attribute (cmesh, itree, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree,
                                sizeof (int), 1);
      }
      t8_cmesh_set_join (cmesh, 0, 1, 1, 0, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
      t8_cmesh_set_partition_range (cmesh, 3, 0, 1);
    }
    else if (rank == 1) {
      t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_QUAD);
      t8_cmesh_set_tree_vertices (cmesh, 1, cells[1].data (), cells[1].size ());
      t8_cmesh_set_tree_vertices (cmesh, 2, cells[2].data (), cells[2].size ());
      t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
      t8_cmesh_set_attribute (cmesh, 2, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
      t8_cmesh_set_join (cmesh, 1, 2, 1, 0, 0);
      t8_cmesh_set_partition_range (cmesh, 3, 2, 2);
    }
  }
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh, 2);
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  t8_forest_t forest
    = t8_forest_new_uniform_with_cmesh_partition (cmesh, t8_scheme_new_default_cxx (), 0, 0, sc_MPI_COMM_WORLD);
  // The following line fails, as the refinement cmesh partitioned required from t8code does not match the provided partition.
  //  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), 0, 1, sc_MPI_COMM_WORLD);
  t8_forest_unref (&forest);
}