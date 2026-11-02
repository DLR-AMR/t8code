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

/** \file t8_gtest_touching_ranks.cxx
 * Tests the ranks touching the part of a rank (\ref t8_forest_touching_ranks.hxx). Two parts touch if and only if
 * two of their leaves touch, so the touching ranks must be the ranks that get ghosts of this rank. They must also be
 * symmetric, have the ranges of the offset arrays of all processes, and be inherited by adapted forests. The forests
 * are uniform, adapted and partitioned, and gathered on one rank, on hypercubes of one or several trees, also
 * periodic ones, in which a neighbor of a cover element may touch it at two opposite boundaries.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_touching_ranks.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_scheme.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include <algorithm>
#include <utility>
#include <vector>

/* The pairs (rank, touching rank) of all ranks */
static std::vector<std::pair<int, int>>
t8_test_allgather_pairs (const std::vector<std::pair<int, int>> &local, sc_MPI_Comm comm)
{
  int size;
  int mpiret = sc_MPI_Comm_size (comm, &size);
  SC_CHECK_MPI (mpiret);
  std::vector<int> flat;
  for (const auto &[a, b] : local) {
    flat.push_back (a);
    flat.push_back (b);
  }
  int count = flat.size ();
  std::vector<int> counts (size), displs (size + 1, 0);
  mpiret = sc_MPI_Allgather (&count, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  for (int rank = 0; rank < size; rank++) {
    displs[rank + 1] = displs[rank] + counts[rank];
  }
  std::vector<int> all_flat (displs[size]);
  mpiret = sc_MPI_Allgatherv (flat.data (), count, sc_MPI_INT, all_flat.data (), counts.data (), displs.data (),
                              sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  std::vector<std::pair<int, int>> all;
  for (size_t i = 0; i < all_flat.size (); i += 2) {
    all.emplace_back (all_flat[i], all_flat[i + 1]);
  }
  return all;
}

/* Check the touching ranks of a forest with a ghost layer against its remote ranks, for symmetry, and against the
 * offset arrays of all processes */
static void
t8_test_check_touching_ranks (const t8_forest_t forest, const int min_dim)
{
  /* The offsets of all processes, which forests only compute on request, as reference */
  if (forest->element_offsets == NULL) {
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == NULL) {
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    t8_forest_partition_create_first_desc (forest);
  }
  ASSERT_NE (forest->touching_ranks, nullptr);
  EXPECT_EQ (forest->touching_ranks->min_dim, min_dim);
  const int rank = forest->mpirank, size = forest->mpisize;

  /* The ranks with ghosts of this rank are the touching ranks */
  std::vector<int> touching, remotes;
  for (const t8_forest_rank_range &range : forest->touching_ranks->ranges) {
    if (range.rank != rank) {
      touching.push_back (range.rank);
    }
  }
  if (forest->ghosts != NULL) {
    int num_remotes;
    const int *remote_ranks = t8_forest_ghost_get_remotes (forest, &num_remotes);
    remotes.assign (remote_ranks, remote_ranks + num_remotes);
    std::sort (remotes.begin (), remotes.end ());
  }
  EXPECT_EQ (touching, remotes) << "touching ranks differ from the ranks with ghosts";

  /* Touching is symmetric */
  std::vector<std::pair<int, int>> pairs;
  for (const int other : touching) {
    pairs.emplace_back (rank, other);
  }
  const std::vector<std::pair<int, int>> all_pairs = t8_test_allgather_pairs (pairs, t8_forest_get_mpicomm (forest));
  for (const auto &[a, b] : all_pairs) {
    if (b == rank) {
      EXPECT_TRUE (std::find (touching.begin (), touching.end (), a) != touching.end ())
        << "rank " << a << " touches this rank, but not vice versa";
    }
  }

  /* The ranges are those of the offset arrays */
  const t8_gloidx_t *offsets = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const t8_gloidx_t *tree_offsets = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
  const t8_linearidx_t *first_descs = (const t8_linearidx_t *) t8_shmem_array_get_array (forest->global_first_desc);
  const auto has_elements = [&] (const int r) { return offsets[r] < offsets[r + 1]; };
  const auto position = [&] (const int r) {
    return t8_forest_position { t8_offset_first (r, tree_offsets), first_descs[r] };
  };
  for (const t8_forest_rank_range &range : forest->touching_ranks->ranges) {
    ASSERT_TRUE (has_elements (range.rank));
    EXPECT_TRUE (range.first == position (range.rank)) << "first position of rank " << range.rank;
    int next = range.rank + 1;
    while (next < size && !has_elements (next)) {
      next++;
    }
    const t8_forest_position next_position
      = next < size ? position (next) : t8_forest_position { t8_forest_get_num_global_trees (forest), 0 };
    EXPECT_TRUE (range.next == next_position) << "next position of rank " << range.rank;
  }
}

/* Refine the elements with (tree + linear id) divisible by 3 down to the level in the user data */
static int
t8_test_touching_adapt (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                        const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                        const t8_scheme_c *scheme, [[maybe_unused]] const int is_family,
                        [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const int max_level = *(const int *) t8_forest_get_user_data (forest);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
  const t8_linearidx_t id = scheme->element_get_linear_id (tree_class, elements[0], level);
  return level < max_level && (gtreeid + id) % 3 == 0;
}

struct forest_touching_ranks: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, bool>>
{
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    eclass = std::get<1> (std::get<0> (GetParam ()));
    const bool periodic = std::get<1> (GetParam ());
    int size;
    const int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
    SC_CHECK_MPI (mpiret);
    if (size == 1) {
      GTEST_SKIP () << "No ghost layer on a single process";
    }
    if (eclass == T8_ECLASS_VERTEX || (periodic && eclass == T8_ECLASS_PYRAMID)) {
      GTEST_SKIP () << "No boundaries, or no periodic hypercube";
    }
    scheme = create_from_scheme_id (scheme_id);
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, periodic);
    const int dim = t8_eclass_to_dimension[eclass];
    uniform_level = dim == 3 ? 1 : dim == 2 ? 2 : 3;
    /* Edges and vertices need the boundaries of the standalone scheme */
    ghost_types.push_back (T8_GHOST_FACES);
    if (scheme_id == 1 && dim >= 2) {
      ghost_types.push_back (T8_GHOST_VERTICES);
      if (dim == 3) {
        ghost_types.push_back (T8_GHOST_EDGES);
      }
    }
  }
  void
  TearDown () override
  {
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
      scheme->unref ();
    }
  }

  /* Check the touching ranks of a forest with a ghost layer of each type, and of an adapted forest that inherits them */
  void
  check_forest (const t8_forest_t forest)
  {
    const int dim = t8_eclass_to_dimension[eclass];
    for (const t8_ghost_type_t ghost_type : ghost_types) {
      const int min_dim = ghost_type == T8_GHOST_FACES ? dim - 1 : ghost_type == T8_GHOST_EDGES ? dim - 2 : 0;
      /* The copy has no touching ranks yet, they are computed for the ghost layer */
      ASSERT_EQ (forest->touching_ranks, nullptr);
      t8_forest_ref (forest);
      t8_forest_t forest_ghost;
      t8_forest_init (&forest_ghost);
      t8_forest_set_copy (forest_ghost, forest);
      t8_forest_set_ghost (forest_ghost, 1, ghost_type);
      t8_forest_commit (forest_ghost);
      t8_test_check_touching_ranks (forest_ghost, min_dim);

      /* Adapting does not change the parts of the ranks */
      int max_level = uniform_level + 2;
      t8_forest_ref (forest_ghost);
      t8_forest_t forest_adapt;
      t8_forest_init (&forest_adapt);
      t8_forest_set_user_data (forest_adapt, &max_level);
      t8_forest_set_adapt (forest_adapt, forest_ghost, t8_test_touching_adapt, 1);
      t8_forest_set_ghost (forest_adapt, 1, ghost_type);
      t8_forest_commit (forest_adapt);
      ASSERT_NE (forest_adapt->touching_ranks, nullptr);
      EXPECT_NE (forest_adapt->touching_ranks, forest_ghost->touching_ranks);
      t8_test_check_touching_ranks (forest_adapt, min_dim);
      t8_forest_unref (&forest_adapt);
      t8_forest_unref (&forest_ghost);
    }
  }

  const t8_scheme *scheme = NULL;
  t8_eclass_t eclass;
  t8_cmesh_t cmesh = NULL;
  int uniform_level;
  std::vector<t8_ghost_type_t> ghost_types;
};

TEST_P (forest_touching_ranks, uniform_adapted_gathered)
{
  /* Level 0 has few elements, so that processes without elements occur */
  for (const int level : { 0, uniform_level }) {
    t8_cmesh_ref (cmesh);
    scheme->ref ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    check_forest (forest);
    t8_forest_unref (&forest);
  }

  /* Adapted and partitioned */
  t8_cmesh_ref (cmesh);
  scheme->ref ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, uniform_level, 0, sc_MPI_COMM_WORLD);
  int max_level = uniform_level + 2;
  t8_forest_t forest_adapt;
  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &max_level);
  t8_forest_set_adapt (forest_adapt, forest, t8_test_touching_adapt, 1);
  t8_forest_set_partition (forest_adapt, NULL, 0);
  t8_forest_commit (forest_adapt);
  check_forest (forest_adapt);

  /* Gathered on the first rank, all other processes have no elements */
  t8_forest_ref (forest_adapt);
  t8_forest_t forest_gather = t8_forest_new_gather (forest_adapt, 0);
  check_forest (forest_gather);
  t8_forest_unref (&forest_gather);
  t8_forest_unref (&forest_adapt);
}

auto print_touching_ranks_param
  = [] (const testing::TestParamInfo<std::tuple<std::tuple<int, t8_eclass_t>, bool>> &info) {
      const int scheme_id = std::get<0> (std::get<0> (info.param));
      const t8_eclass_t eclass = std::get<1> (std::get<0> (info.param));
      return std::string (t8_scheme_to_string[scheme_id]) + "_" + t8_eclass_to_string[eclass]
             + (std::get<1> (info.param) ? "_periodic" : "");
    };

INSTANTIATE_TEST_SUITE_P (t8_gtest_touching_ranks, forest_touching_ranks,
                          ::testing::Combine (AllSchemes, ::testing::Bool ()), print_touching_ranks_param);
