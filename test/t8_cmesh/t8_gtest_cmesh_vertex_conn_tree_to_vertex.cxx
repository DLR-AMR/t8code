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
#include <test/t8_gtest_macros.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity.hxx>
#include <t8_cmesh/t8_cmesh_vertex_conn_tree_to_vertex.hxx>
#include <t8_cmesh/t8_cmesh_vertex_conn_vertex_to_tree.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <test/t8_cmesh_generator/t8_cmesh_example_sets.hxx>

/* TODO: write test case without existing cmesh to test before attribute bug is fixed */
class cmesh_vertex_conn_ttv_with_core_classes: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    const t8_cmesh_t committed_cmesh = GetParam ()->cmesh_create ();
    T8_ASSERT (t8_cmesh_is_committed (committed_cmesh));
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_derive (cmesh, committed_cmesh);
    t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
    t8_cmesh_set_partition_uniform (cmesh, 0, scheme);
    const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (committed_cmesh);

    t8_debugf ("Starting test with cmesh of dim %i and %li global, %i local trees.\n", cmesh->dimension,
               t8_cmesh_get_num_trees (committed_cmesh), num_local_trees);
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

      const t8_eclass_t tree_class = t8_cmesh_get_tree_class (committed_cmesh, itree);
      /* Allocate space for entries. */
      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];
      t8_gloidx_t *global_indices = T8_ALLOC (t8_gloidx_t, num_tree_vertices);
      /* Fill with 0, 1, 2, 3, 4 ... */
      const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
      for (t8_locidx_t ientry = start_index; ientry < num_tree_vertices; ++ientry) {
        global_indices[ientry] = ientry;
      }
      ttv.set_global_vertex_ids_of_tree_vertices (cmesh, itree, global_indices, num_tree_vertices);
      /* It is save to free the entries after commit, since the value got copied. */
      T8_FREE (global_indices);
    }

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;

  t8_cmesh_vertex_conn_tree_to_vertex ttv;
};

/** Check attribute values of cmeshes against reference values. */
TEST_P (cmesh_vertex_conn_ttv_with_core_classes, DISABLED_get_global)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Get all vertices */
    const t8_gloidx_t *global_vertices = ttv.get_global_vertices (cmesh, itree, num_tree_vertices);

    const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      /* Get the stored global vertex id */
      const t8_gloidx_t global_vertex = ttv.get_global_vertex (cmesh, itree, ivertex, num_tree_vertices);
      /* Check value */
      EXPECT_EQ (global_vertex, start_index + ivertex);
      EXPECT_EQ (global_vertices[ivertex], start_index + ivertex);
    }
  }
}

/* TODO: Remove the tests belonging to cmesh_vertex_conn_ttv_temp
 *       as soon as we can enable the tests cmesh_vertex_conn_ttv.
 *       That is as soon as we can add attributes to cmeshes while deriving. */

#define VTT_TEST_MAX_NUM_TREES 100

class cmesh_vertex_conn_ttv_with_core_classes_temp:
  public testing::TestWithParam<std::tuple<t8_gloidx_t, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    num_trees = std::get<0> (GetParam ());
    const t8_eclass_t tree_class = std::get<1> (GetParam ());
    t8_cmesh_init (&cmesh);

    t8_debugf ("Testing cmesh with %li trees of class %s\n", num_trees, t8_eclass_to_string[tree_class]);

    for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
      /* Set this tree's class. */
      t8_cmesh_set_tree_class (cmesh, itree, tree_class);
      /* Allocate space for entries. */
      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];
      t8_gloidx_t *global_indices = T8_ALLOC (t8_gloidx_t, num_tree_vertices);
      /* Fill with i, i+1, i+2, i+3, ... */
      const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
      for (t8_locidx_t ientry = 0; ientry < num_tree_vertices; ++ientry) {
        global_indices[ientry] = start_index + ientry;
      }
      ttv.set_global_vertex_ids_of_tree_vertices (cmesh, itree, global_indices, num_tree_vertices);
      /* It is save to free the entries after commit, since the value got copied. */
      T8_FREE (global_indices);
    }

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_gloidx_t num_trees;

  t8_cmesh_vertex_conn_tree_to_vertex ttv;
};

/** Check for correct ttv entries. */
#if 0
// Reactive this line when we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv, DISABLED_get_global)
#else
// Delete this line and the cmesh_vertex_conn_ttv_temp class wehen we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv_with_core_classes_temp, get_global)
#endif
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Get all vertices */
    const t8_gloidx_t *global_vertices = ttv.get_global_vertices (cmesh, itree, num_tree_vertices);

    const t8_gloidx_t global_tree_id = t8_cmesh_get_global_id (cmesh, itree);
    const t8_gloidx_t start_index = global_tree_id * T8_ECLASS_MAX_CORNERS;
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      /* Get the stored global vertex id */
      const t8_gloidx_t global_vertex = ttv.get_global_vertex (cmesh, itree, ivertex, num_tree_vertices);
      /* Check value */
      EXPECT_EQ (global_vertex, start_index + ivertex);
      EXPECT_EQ (global_vertices[ivertex], start_index + ivertex);
    }
  }
}

#if 0
// Reactive this line when we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv, DISABLED_get_global)
#else
// Delete this line and the cmesh_vertex_conn_ttv_temp class wehen we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv_with_core_classes_temp, convert_to_vtt)
#endif
{
  /* create a vertex_to_tree list from ttv */
  t8_cmesh_vertex_conn_vertex_to_tree vtt;
  vtt.build_from_ttv(cmesh, ttv);
  /* Since global tree i is mapped to vertices:
   *  i*T8_ECLASS_MAX_CORNERS, i*T8_ECLASS_MAX_CORNERS + 1, ... 
   *  and this mapping is unique, we know that the list for vertex j
   *  must contain 
   *  global tree j / T8_ECLASS_MAX_CORNERS
   *  with local vertex j % T8_ECLASS_MAX_CORNERS */
  ASSERT_TRUE (vtt.is_committed ());

  /* Iterate over all entries in vtt.
   * Each entry corresponds to a global vertex id and
   * gives its list of tree indices and vertices. */
  for (auto &[global_vertex, tree_vertex_list] : vtt) {
    /* Iterate over the list of tree indices and vertices of this global vertex. */
    for (auto &[tree_index, tree_vertex] : tree_vertex_list) {
      const t8_locidx_t expected_tree = global_vertex / T8_ECLASS_MAX_CORNERS;
      const int expected_tree_vertex = global_vertex % T8_ECLASS_MAX_CORNERS;
      EXPECT_EQ (tree_index, expected_tree);
      EXPECT_EQ (tree_vertex, expected_tree_vertex);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_tree_to_vertex, cmesh_vertex_conn_ttv_with_core_classes, AllCmeshsParam,
                          pretty_print_base_example);

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_tree_to_vertex, cmesh_vertex_conn_ttv_with_core_classes_temp,
                          testing::Combine (testing::Values (1, VTT_TEST_MAX_NUM_TREES + 1), AllEclasses));

class cmesh_vertex_conn_ttv_with_cmesh_functions: public testing::TestWithParam<cmesh_example_base *> {
 protected:
  void
  SetUp () override
  {
    const t8_cmesh_t committed_cmesh = GetParam ()->cmesh_create ();
    T8_ASSERT (t8_cmesh_is_committed (committed_cmesh));
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_derive (cmesh, committed_cmesh);
    t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();
    t8_cmesh_set_partition_uniform (cmesh, 0, scheme);
    const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (committed_cmesh);

    t8_debugf ("Starting test with cmesh of dim %i and %li global, %i local trees.\n", cmesh->dimension,
               t8_cmesh_get_num_trees (committed_cmesh), num_local_trees);
    for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

      const t8_eclass_t tree_class = t8_cmesh_get_tree_class (committed_cmesh, itree);
      /* Allocate space for entries. */
      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];
      t8_gloidx_t *global_indices = T8_ALLOC (t8_gloidx_t, num_tree_vertices);
      /* Fill with 0, 1, 2, 3, 4 ... */
      const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
      for (t8_locidx_t ientry = start_index; ientry < num_tree_vertices; ++ientry) {
        global_indices[ientry] = ientry;
      }
      t8_cmesh_set_global_vertices_of_tree (cmesh, itree, global_indices, num_tree_vertices);
      /* It is save to free the entries after commit, since the value got copied. */
      T8_FREE (global_indices);
    }

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;
};

/** Check attribute values of cmeshes against reference values. */
TEST_P (cmesh_vertex_conn_ttv_with_cmesh_functions, DISABLED_get_global)
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Get all vertices */
    const t8_gloidx_t *global_vertices = t8_cmesh_get_global_vertices_of_tree (cmesh, itree, num_tree_vertices);

    const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      /* Get the stored global vertex id */
      const t8_gloidx_t global_vertex = t8_cmesh_get_global_vertex_of_tree (cmesh, itree, ivertex, num_tree_vertices);
      /* Check value */
      EXPECT_EQ (global_vertex, start_index + ivertex);
      EXPECT_EQ (global_vertices[ivertex], start_index + ivertex);
    }
  }
}

/* TODO: Remove the tests belonging to cmesh_vertex_conn_ttv_temp
 *       as soon as we can enable the tests cmesh_vertex_conn_ttv.
 *       That is as soon as we can add attributes to cmeshes while deriving. */

#define VTT_TEST_MAX_NUM_TREES 100

class cmesh_vertex_conn_ttv_with_cmesh_functions_temp:
  public testing::TestWithParam<std::tuple<t8_gloidx_t, t8_eclass_t>> {
 protected:
  void
  SetUp () override
  {
    num_trees = std::get<0> (GetParam ());
    const t8_eclass_t tree_class = std::get<1> (GetParam ());
    t8_cmesh_init (&cmesh);

    t8_debugf ("Testing cmesh with %li trees of class %s\n", num_trees, t8_eclass_to_string[tree_class]);

    for (t8_locidx_t itree = 0; itree < num_trees; ++itree) {
      /* Set this tree's class. */
      t8_cmesh_set_tree_class (cmesh, itree, tree_class);
      /* Allocate space for entries. */
      const int num_tree_vertices = t8_eclass_num_vertices[tree_class];
      t8_gloidx_t *global_indices = T8_ALLOC (t8_gloidx_t, num_tree_vertices);
      /* Fill with i, i+1, i+2, i+3, ... */
      const t8_gloidx_t start_index = itree * T8_ECLASS_MAX_CORNERS;
      for (t8_locidx_t ientry = 0; ientry < num_tree_vertices; ++ientry) {
        global_indices[ientry] = start_index + ientry;
      }
      t8_cmesh_set_global_vertices_of_tree (cmesh, itree, global_indices, num_tree_vertices);
      /* It is save to free the entries after commit, since the value got copied. */
      T8_FREE (global_indices);
    }

    /* Commit the cmesh */
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_cmesh_t cmesh;
  t8_gloidx_t num_trees;
};

/** Check for correct ttv entries. */
#if 0
// Reactive this line when we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv, DISABLED_get_global)
#else
// Delete this line and the cmesh_vertex_conn_ttv_temp class wehen we enable the tests with derived attributes
TEST_P (cmesh_vertex_conn_ttv_with_cmesh_functions_temp, get_global)
#endif
{
  const t8_locidx_t num_local_trees = t8_cmesh_get_num_local_trees (cmesh);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {

    const t8_eclass_t tree_class = t8_cmesh_get_tree_class (cmesh, itree);
    const int num_tree_vertices = t8_eclass_num_vertices[tree_class];

    /* Get all vertices */
    const t8_gloidx_t *global_vertices = t8_cmesh_get_global_vertices_of_tree (cmesh, itree, num_tree_vertices);

    const t8_gloidx_t global_tree_id = t8_cmesh_get_global_id (cmesh, itree);
    const t8_gloidx_t start_index = global_tree_id * T8_ECLASS_MAX_CORNERS;
    for (int ivertex = 0; ivertex < num_tree_vertices; ++ivertex) {
      /* Get the stored global vertex id */
      const t8_gloidx_t global_vertex = t8_cmesh_get_global_vertex_of_tree (cmesh, itree, ivertex, num_tree_vertices);
      /* Check value */
      EXPECT_EQ (global_vertex, start_index + ivertex);
      EXPECT_EQ (global_vertices[ivertex], start_index + ivertex);
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_tree_to_vertex, cmesh_vertex_conn_ttv_with_cmesh_functions,
                          AllCmeshsParam, pretty_print_base_example);

INSTANTIATE_TEST_SUITE_P (t8_gtest_cmesh_vertex_tree_to_vertex, cmesh_vertex_conn_ttv_with_cmesh_functions_temp,
                          testing::Combine (testing::Values (1, VTT_TEST_MAX_NUM_TREES + 1), AllEclasses));