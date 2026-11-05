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

/** \file t8_gtest_shmem.cxx
 * Tests of the shared memory array \ref t8_shmem_array and of its C interface.
 * In debugging mode, every test runs on several simulated node layouts (see T8_SHMEM_SIMULATE_NODES in
 * t8_shmem.cxx), so that the exchange between nodes is tested on a single machine as well.
 */

#include <gtest/gtest.h>
#include <t8_data/t8_shmem.hxx>
#include <test/t8_gtest_memory_macros.hxx>

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

static_assert (std::ranges::contiguous_range<t8_shmem_array<t8_gloidx_t>>);
static_assert (std::ranges::sized_range<t8_shmem_array<t8_gloidx_t>>);
static_assert (std::is_nothrow_move_constructible_v<t8_shmem_array<t8_gloidx_t>>);
static_assert (!std::is_copy_constructible_v<t8_shmem_array<t8_gloidx_t>>);

/** Parameter: the communicator to duplicate and the simulated node layout ("" for the real layout). */
struct shmem: public testing::TestWithParam<std::tuple<sc_MPI_Comm, const char *>>
{
 protected:
  void
  SetUp () override
  {
    layout = std::get<1> (GetParam ());
    int mpiret = sc_MPI_Comm_dup (std::get<0> (GetParam ()), &comm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    /* The layout is read when the topology of the fresh communicator is computed. */
    if (!layout.empty ()) {
      setenv ("T8_SHMEM_SIMULATE_NODES", layout.c_str (), 1);
    }
    shared_size = t8_shmem_init (comm);
    unsetenv ("T8_SHMEM_SIMULATE_NODES");
  }

  void
  TearDown () override
  {
    t8_shmem_finalize (comm);
    const int mpiret = sc_MPI_Comm_free (&comm);
    SC_CHECK_MPI (mpiret);
  }

  /** Check that all processes of the communicator see the same content. */
  template <typename TValue>
  void
  check_identical (const t8_shmem_array<TValue> &array)
  {
    const std::vector<TValue> local (array.begin (), array.end ());
    std::vector<TValue> root (local);
    const int mpiret
      = sc_MPI_Bcast (root.data (), static_cast<int> (root.size () * sizeof (TValue)), sc_MPI_BYTE, 0, comm);
    SC_CHECK_MPI (mpiret);
    EXPECT_EQ (local, root);
  }

  sc_MPI_Comm comm;   /**< A duplicate of the parameter communicator with its own topology. */
  int mpirank;        /**< The rank in \a comm. */
  int mpisize;        /**< The size of \a comm. */
  int shared_size;    /**< The return value of t8_shmem_init. */
  std::string layout; /**< The simulated node layout. */
};

TEST_P (shmem, test_init_and_mode)
{
  EXPECT_GE (shared_size, 1);
  t8_shmem_array<t8_gloidx_t> array (comm, mpisize + 1);
  EXPECT_EQ (array.get_comm (), comm);
  EXPECT_EQ (array.size (), (size_t) mpisize + 1);
  /* A node with a single process still uses a window if other nodes share memory. */
  if (shared_size > 1) {
    EXPECT_TRUE (array.is_node_shared ());
  }
#if T8_ENABLE_DEBUG && T8_ENABLE_MPI && defined(SC_ENABLE_MPIWINSHARED) && defined(SC_ENABLE_MPICOMMSHARED)
  /* Replicated, non-contiguous and single-process nodes cannot share memory. */
  if (layout == "replicated" || layout == "block:1" || (layout == "cyclic:2" && mpisize > 2)) {
    EXPECT_FALSE (array.is_node_shared ());
  }
  if (layout == "block:2" && mpisize >= 2) {
    EXPECT_TRUE (array.is_node_shared ());
  }
#endif
}

TEST_P (shmem, test_allgather)
{
  t8_shmem_array<t8_gloidx_t> array (comm, mpisize + 1);
  array.allgather (10 * mpirank, [this] (std::span<t8_gloidx_t> data) { data[mpisize] = -1; });
  for (int iproc = 0; iproc < mpisize; ++iproc) {
    EXPECT_EQ (array[iproc], 10 * iproc);
  }
  EXPECT_EQ (array.back (), -1);
  check_identical (array);

  /* Refill after every process has read the array. */
  const t8_gloidx_t sum = std::accumulate (array.begin (), array.end (), t8_gloidx_t { 0 });
  EXPECT_EQ (sum, 5 * mpisize * (t8_gloidx_t) (mpisize - 1) - 1);
  array.allgather (sum + mpirank);
  for (int iproc = 0; iproc < mpisize; ++iproc) {
    EXPECT_EQ (array[iproc], sum + iproc);
  }
  EXPECT_EQ (array.back (), -1);

  /* Other element types */
  t8_shmem_array<t8_linearidx_t> linear_ids (comm, mpisize);
  linear_ids.allgather ((t8_linearidx_t) mpirank << 40);
  for (int iproc = 0; iproc < mpisize; ++iproc) {
    EXPECT_EQ (linear_ids[iproc], (t8_linearidx_t) iproc << 40);
  }
}

TEST_P (shmem, test_allgatherv)
{
  /* Rank r sends r % 3 values, so some processes send nothing. */
  std::vector<t8_gloidx_t> expected;
  for (int iproc = 0; iproc < mpisize; ++iproc) {
    for (int ivalue = 0; ivalue < iproc % 3; ++ivalue) {
      expected.push_back (100 * iproc + ivalue);
    }
  }
  std::vector<t8_gloidx_t> local;
  for (int ivalue = 0; ivalue < mpirank % 3; ++ivalue) {
    local.push_back (100 * mpirank + ivalue);
  }

  /* Two extra elements behind the gathered data must stay untouched. */
  t8_shmem_array<t8_gloidx_t> array (comm, expected.size () + 2);
  array.write ([] (std::span<t8_gloidx_t> data) { std::ranges::fill (data, 77); });
  array.allgatherv (local);
  EXPECT_TRUE (std::equal (expected.begin (), expected.end (), array.begin ()));
  EXPECT_EQ (array[expected.size ()], 77);
  EXPECT_EQ (array.back (), 77);
  check_identical (array);

  /* Refill with a finalize that patches the first entry. */
  array.allgatherv (local, [] (std::span<t8_gloidx_t> data) { data[0] = -5; });
  EXPECT_EQ (array.front (), -5);
  EXPECT_TRUE (std::equal (expected.begin () + std::min<size_t> (1, expected.size ()), expected.end (),
                           array.begin () + std::min<size_t> (1, expected.size ())));

  /* Bytes with an empty total */
  t8_shmem_array<char> empty (comm, 0);
  empty.allgatherv (std::span<const char> ());
  EXPECT_TRUE (empty.empty ());
  EXPECT_EQ (empty.begin (), empty.end ());
}

TEST_P (shmem, test_prefix_sum)
{
  t8_shmem_array<t8_gloidx_t> offsets (comm, mpisize + 1);
  offsets.prefix_sum (mpirank + 1);
  t8_gloidx_t expected = 0;
  for (int iproc = 0; iproc <= mpisize; ++iproc) {
    EXPECT_EQ (offsets[iproc], expected);
    expected += iproc + 1;
  }
  check_identical (offsets);

  /* Owner search with std algorithms: the owner p of a value v satisfies offsets[p] <= v < offsets[p + 1]. */
  for (t8_gloidx_t value = 0; value < offsets.back (); ++value) {
    const int owner = std::ranges::upper_bound (offsets, value) - offsets.begin () - 1;
    EXPECT_LE (offsets[owner], value);
    EXPECT_LT (value, offsets[owner + 1]);
  }
}

TEST_P (shmem, test_write_clone_compare_move)
{
  t8_shmem_array<t8_gloidx_t> array (comm, 2 * mpisize);
  array.write ([] (std::span<t8_gloidx_t> data) { std::iota (data.begin (), data.end (), 0); });
  /* Writing a published array */
  array.write ([] (std::span<t8_gloidx_t> data) {
    for (auto &value : data) {
      value *= 2;
    }
  });
  for (size_t index = 0; index < array.size (); ++index) {
    EXPECT_EQ (array[index], 2 * (t8_gloidx_t) index);
  }

  t8_shmem_array<t8_gloidx_t> copy = array.clone ();
  EXPECT_TRUE (copy == array);
  EXPECT_NE (copy.data (), array.data ());
  copy.write ([] (std::span<t8_gloidx_t> data) { data[0] = -1; });
  EXPECT_FALSE (copy == array);

  t8_shmem_array<t8_gloidx_t> moved (std::move (copy));
  EXPECT_EQ (moved.front (), -1);
  EXPECT_TRUE (copy.empty ());
  moved = array.clone ();
  EXPECT_TRUE (moved == array);

  const std::span<const t8_gloidx_t> view = array;
  EXPECT_EQ (view.size (), array.size ());
  EXPECT_EQ (view.data (), array.data ());
}

TEST_P (shmem, test_c_interface)
{
  t8_shmem_array_t array;
  t8_shmem_array_init (&array, sizeof (t8_gloidx_t), mpisize + 1, comm);
  EXPECT_EQ (t8_shmem_array_get_elem_count (array), (size_t) mpisize + 1);
  EXPECT_EQ (t8_shmem_array_get_elem_size (array), sizeof (t8_gloidx_t));
  EXPECT_EQ (t8_shmem_array_get_comm (array), comm);

  const t8_gloidx_t value = mpirank;
  t8_shmem_array_allgather (&value, 1, T8_MPI_GLOIDX, array, 1, T8_MPI_GLOIDX);
  if (t8_shmem_array_start_writing (array)) {
    t8_shmem_array_set_gloidx (array, mpisize, mpisize);
    *(t8_gloidx_t *) t8_shmem_array_index_for_writing (array, 0) = 0;
  }
  t8_shmem_array_end_writing (array);
  for (int iproc = 0; iproc <= mpisize; ++iproc) {
    EXPECT_EQ (t8_shmem_array_get_gloidx (array, iproc), iproc);
    EXPECT_EQ (*(const t8_gloidx_t *) t8_shmem_array_index (array, iproc), iproc);
    EXPECT_EQ (t8_shmem_array_get_gloidx_array (array)[iproc], iproc);
  }

  t8_shmem_array_t copy;
  t8_shmem_array_init (&copy, sizeof (t8_gloidx_t), mpisize + 1, comm);
  t8_shmem_array_copy (copy, array);
  EXPECT_TRUE (t8_shmem_array_is_equal (copy, array));
  EXPECT_FALSE (t8_shmem_array_is_equal (copy, nullptr));

  /* The prefix sum of ones equals the identity. */
  const t8_gloidx_t one = 1;
  t8_shmem_array_prefix (&one, copy, 1, T8_MPI_GLOIDX, sc_MPI_SUM, comm);
  EXPECT_TRUE (t8_shmem_array_is_equal (copy, array));

  /* Every process sends one value, rank p sends p. */
  t8_shmem_array_allgatherv ((void *) &value, 1, T8_MPI_GLOIDX, copy, T8_MPI_GLOIDX, comm);
  for (int iproc = 0; iproc < mpisize; ++iproc) {
    EXPECT_EQ (t8_shmem_array_get_gloidx (copy, iproc), iproc);
  }

  auto compare = [] (t8_shmem_array_t shmem_array, const int index, const t8_gloidx_t search) {
    const t8_gloidx_t entry = t8_shmem_array_get_gloidx (shmem_array, index);
    return search == entry ? 0 : (search < entry ? -1 : 1);
  };
  for (int iproc = 0; iproc <= mpisize; ++iproc) {
    EXPECT_EQ (t8_shmem_array_binary_search (array, iproc, mpisize, compare), iproc);
  }
  EXPECT_EQ (t8_shmem_array_binary_search (array, mpisize + 1, mpisize, compare), -1);

  t8_shmem_array_destroy (&copy);
  t8_shmem_array_destroy (&array);
  EXPECT_EQ (array, nullptr);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_shmem, shmem,
                          testing::Combine (testing::Values (sc_MPI_COMM_WORLD, sc_MPI_COMM_SELF),
                                            testing::Values ("", "replicated", "block:1", "block:2", "block:3",
                                                             "cyclic:2")));
