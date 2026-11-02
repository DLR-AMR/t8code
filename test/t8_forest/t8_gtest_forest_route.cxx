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

/** \file t8_gtest_forest_route.cxx
 * Tests the partition of a forest as seen by one process (\ref t8_forest_splitters.hxx) against the offset arrays of
 * all processes, and the routing of items along it: every item must arrive exactly at the processes that own a
 * position between its first and last position (or a global element index in its range), and with \a unique only at
 * its single owner. The forests are uniform, adapted and partitioned, and gathered on one rank, so that processes
 * without elements occur. The reference owners are computed from the offset arrays of all processes.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_offset.h>
#include <t8_data/t8_shmem.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_splitters.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_scheme.hxx>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include <algorithm>
#include <random>
#include <vector>

/* An element of a tree, routed to the owners of its descendants */
struct t8_test_route_item
{
  t8_forest_position first; /* The position of the first descendant */
  t8_forest_position last;  /* The position of the last descendant */
  t8_gloidx_t gtreeid;      /* The tree of the element */
  t8_linearidx_t id;        /* The linear id of the element at its level */
  int level;                /* The level of the element */
  int source;               /* The rank that created the item */

  friend auto
  operator<=> (const t8_test_route_item &, const t8_test_route_item &)
    = default;
};

/* A range of global element indices, routed to the owners of its elements */
struct t8_test_index_item
{
  t8_gloidx_t first; /* The first global element index */
  t8_gloidx_t last;  /* The last global element index */
  int source;        /* The rank that created the item */
  int serial;        /* The number of the item on its source */

  friend auto
  operator<=> (const t8_test_index_item &, const t8_test_index_item &)
    = default;
};

/* The items of all processes, in the order of the ranks */
template <typename T>
static std::vector<T>
t8_test_allgather (const std::vector<T> &local, sc_MPI_Comm comm)
{
  int size;
  int mpiret = sc_MPI_Comm_size (comm, &size);
  SC_CHECK_MPI (mpiret);
  int num_bytes = local.size () * sizeof (T);
  std::vector<int> counts (size), displs (size + 1, 0);
  mpiret = sc_MPI_Allgather (&num_bytes, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  for (int rank = 0; rank < size; rank++) {
    displs[rank + 1] = displs[rank] + counts[rank];
  }
  std::vector<T> all (displs[size] / sizeof (T));
  mpiret = sc_MPI_Allgatherv (const_cast<T *> (local.data ()), num_bytes, sc_MPI_BYTE, all.data (), counts.data (),
                              displs.data (), sc_MPI_BYTE, comm);
  SC_CHECK_MPI (mpiret);
  return all;
}

/* The owners of the positions [first, last] according to the offset arrays of all processes: each process with
 * elements owns the positions from its first element (from the beginning if it is the first process with elements)
 * to the first element of the next process with elements (to the end if there is none). */
static std::vector<int>
t8_test_position_owners (const t8_forest_t forest, const t8_forest_position &first, const t8_forest_position &last)
{
  const t8_gloidx_t *offsets = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  const t8_gloidx_t *tree_offsets = t8_shmem_array_get_gloidx_array (forest->tree_offsets);
  const t8_linearidx_t *first_descs = (const t8_linearidx_t *) t8_shmem_array_get_array (forest->global_first_desc);
  std::vector<int> ranks;
  for (int rank = 0; rank < forest->mpisize; rank++) {
    if (offsets[rank] < offsets[rank + 1]) {
      ranks.push_back (rank);
    }
  }
  const auto position = [&] (const int rank) {
    return t8_forest_position { t8_offset_first (rank, tree_offsets), first_descs[rank] };
  };
  std::vector<int> owners;
  for (size_t irank = 0; irank < ranks.size (); irank++) {
    const bool starts_before_last = irank == 0 || position (ranks[irank]) <= last;
    const bool ends_after_first = irank + 1 == ranks.size () || first < position (ranks[irank + 1]);
    if (starts_before_last && ends_after_first) {
      owners.push_back (ranks[irank]);
    }
  }
  return owners;
}

/* The owners of the global element indices [first, last] according to the element offsets of all processes */
static std::vector<int>
t8_test_index_owners (const t8_forest_t forest, const t8_gloidx_t first, const t8_gloidx_t last)
{
  const t8_gloidx_t *offsets = t8_shmem_array_get_gloidx_array (forest->element_offsets);
  std::vector<int> owners;
  for (int rank = 0; rank < forest->mpisize; rank++) {
    if (offsets[rank] < offsets[rank + 1] && offsets[rank] <= last && first < offsets[rank + 1]) {
      owners.push_back (rank);
    }
  }
  return owners;
}

/* Whether this rank should receive an item with the given owners */
static bool
t8_test_receives (const std::vector<int> &owners, const int rank, const bool unique)
{
  return unique ? owners.size () == 1 && owners[0] == rank
                : std::find (owners.begin (), owners.end (), rank) != owners.end ();
}

/* The offsets of all processes, which forests only compute on request, as reference */
static void
t8_test_create_offsets (t8_forest_t forest)
{
  if (forest->element_offsets == NULL) {
    t8_forest_partition_create_offsets (forest);
  }
  if (forest->tree_offsets == NULL) {
    t8_forest_partition_create_tree_offsets (forest);
  }
  if (forest->global_first_desc == NULL) {
    t8_forest_partition_create_first_desc (forest);
  }
}

/* Check the splitters of a forest and route elements of all levels up to max_item_level and ranges of global element
 * indices along them */
static void
t8_test_splitters_and_route (const t8_forest_t forest, const int max_item_level)
{
  t8_test_create_offsets (forest);
  ASSERT_TRUE (t8_forest_splitters_match_offsets (forest));
  const sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  const int rank = forest->mpirank, size = forest->mpisize;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);

  /* All elements up to max_item_level of all trees, distributed round robin over the ranks */
  std::vector<t8_test_route_item> items;
  long counter = 0;
  for (t8_gloidx_t gtreeid = 0; gtreeid < t8_forest_get_num_global_trees (forest); gtreeid++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, gtreeid);
    t8_element_t *element;
    scheme->element_new (eclass, 1, &element);
    for (int level = 0; level <= max_item_level; level++) {
      const t8_gloidx_t num_elements = scheme->count_leaves_from_root (eclass, level);
      for (t8_gloidx_t id = 0; id < num_elements; id++) {
        if (counter++ % size != rank) {
          continue;
        }
        t8_test_route_item item;
        scheme->element_set_linear_id (eclass, element, level, id);
        t8_forest_element_positions (forest, gtreeid, eclass, element, &item.first, &item.last);
        item.gtreeid = gtreeid;
        item.id = id;
        item.level = level;
        item.source = rank;
        items.push_back (item);
      }
    }
    scheme->element_destroy (eclass, 1, &element);
  }
  const std::vector<t8_test_route_item> all_items = t8_test_allgather (items, comm);

  for (const bool unique : { false, true }) {
    const auto halves = [unique] (const t8_test_route_item &item, const t8_forest_bisection &level) {
      return t8_forest_route_halves_positions (level, item.first, item.last, unique);
    };
    std::vector<t8_test_route_item> expected;
    for (const t8_test_route_item &item : all_items) {
      if (t8_test_receives (t8_test_position_owners (forest, item.first, item.last), rank, unique)) {
        expected.push_back (item);
      }
    }
    std::sort (expected.begin (), expected.end ());

    std::vector<t8_test_route_item> routed = items;
    t8_forest_route (forest, routed, halves);
    std::sort (routed.begin (), routed.end ());
    EXPECT_TRUE (routed == expected) << "routed " << routed.size () << " items, expected " << expected.size ()
                                     << (unique ? " (unique)" : "");

    /* Every rank routes all items, the merge removes the copies of the other ranks */
    const auto merge = [] (std::vector<t8_test_route_item> &merged) {
      std::sort (merged.begin (), merged.end ());
      merged.erase (std::unique (merged.begin (), merged.end ()), merged.end ());
    };
    routed = all_items;
    t8_forest_route (forest, routed, halves, merge);
    /* Without other ranks, there is no round and no merge */
    merge (routed);
    EXPECT_TRUE (routed == expected) << "merged " << routed.size () << " items, expected " << expected.size ()
                                     << (unique ? " (unique)" : "");
  }

  /* Ranges of global element indices of different lengths */
  const t8_gloidx_t num_elements = t8_forest_get_global_num_leaf_elements (forest);
  if (num_elements == 0) {
    return;
  }
  std::mt19937_64 random (rank);
  std::vector<t8_test_index_item> index_items;
  for (int serial = 0; serial < 30; serial++) {
    const t8_gloidx_t max_length = serial % 3 == 0 ? 1 : serial % 3 == 1 ? num_elements / 8 + 1 : num_elements;
    const t8_gloidx_t first = random () % num_elements;
    const t8_gloidx_t last = std::min (num_elements - 1, first + (t8_gloidx_t) (random () % max_length));
    index_items.push_back ({ first, last, rank, serial });
  }
  const std::vector<t8_test_index_item> all_index_items = t8_test_allgather (index_items, comm);
  for (const bool unique : { false, true }) {
    std::vector<t8_test_index_item> expected;
    for (const t8_test_index_item &item : all_index_items) {
      if (t8_test_receives (t8_test_index_owners (forest, item.first, item.last), rank, unique)) {
        expected.push_back (item);
      }
    }
    std::sort (expected.begin (), expected.end ());
    std::vector<t8_test_index_item> routed = index_items;
    t8_forest_route (forest, routed, [unique] (const t8_test_index_item &item, const t8_forest_bisection &level) {
      return t8_forest_route_halves_elements (level, item.first, item.last, unique);
    });
    std::sort (routed.begin (), routed.end ());
    EXPECT_TRUE (routed == expected) << "routed " << routed.size () << " index ranges, expected " << expected.size ()
                                     << (unique ? " (unique)" : "");
  }
}

/* Refine the elements with (tree + linear id) divisible by 3 down to the level in the user data */
static int
t8_test_route_adapt (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                     const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                     [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                     t8_element_t *elements[])
{
  const int max_level = *(const int *) t8_forest_get_user_data (forest);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
  const t8_linearidx_t id = scheme->element_get_linear_id (tree_class, elements[0], level);
  return level < max_level && (gtreeid + id) % 3 == 0;
}

struct forest_route: public testing::TestWithParam<std::tuple<int, t8_eclass_t>>
{
 protected:
  void
  SetUp () override
  {
    scheme = create_from_scheme_id (std::get<0> (GetParam ()));
    eclass = std::get<1> (GetParam ());
    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    const int dim = t8_eclass_to_dimension[eclass];
    uniform_level = dim == 3 ? 1 : dim == 2 ? 2 : 3;
  }
  void
  TearDown () override
  {
    t8_cmesh_unref (&cmesh);
    scheme->unref ();
  }
  const t8_scheme *scheme;
  t8_eclass_t eclass;
  t8_cmesh_t cmesh;
  int uniform_level;
};

TEST_P (forest_route, uniform_adapted_gathered)
{
  /* Level 0 has few elements, so that processes without elements occur */
  for (const int level : { 0, uniform_level }) {
    t8_cmesh_ref (cmesh);
    scheme->ref ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
    t8_test_splitters_and_route (forest, uniform_level + 1);
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
  t8_forest_set_adapt (forest_adapt, forest, t8_test_route_adapt, 1);
  t8_forest_set_partition (forest_adapt, NULL, 0);
  t8_forest_commit (forest_adapt);
  t8_test_splitters_and_route (forest_adapt, uniform_level + 1);

  /* Gathered on the first and on the last rank, all other processes have no elements */
  int size;
  const int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
  SC_CHECK_MPI (mpiret);
  for (const int gather_rank : { 0, size - 1 }) {
    t8_forest_ref (forest_adapt);
    t8_forest_t forest_gather = t8_forest_new_gather (forest_adapt, gather_rank);
    t8_test_splitters_and_route (forest_gather, uniform_level + 1);
    t8_forest_unref (&forest_gather);
  }
  t8_forest_unref (&forest_adapt);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_forest_route, forest_route, AllSchemes, print_all_schemes);
