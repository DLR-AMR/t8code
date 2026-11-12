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

/**
 * \file t8_gtest_data_exchange.cxx
 *
 * Tests for the repartitioning of distributed arrays: the partition pattern, the byte exchange and the
 * repartitioning of vector handlers. Every process computes the same old and new distributions without
 * communication.
 */

#include <gtest/gtest.h>
#include <t8_data/t8_data_exchange.hxx>
#include <t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_enlarged_stdtypes.hxx>
#include <t8_data/t8_pseudo_trees.hxx>
#include <t8_data/t8_vector_handler.hxx>
#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

/** Distributions used to test repartitioning. */
enum t8_test_partition_shape {
  T8_TEST_PARTITION_SHIFT,    /**< Every process boundary moves by three items. */
  T8_TEST_PARTITION_GATHER,   /**< All items move to process 0. */
  T8_TEST_PARTITION_SCATTER,  /**< All items move from process 0 to all processes. */
  T8_TEST_PARTITION_IDENTITY, /**< No item moves. */
  T8_TEST_PARTITION_REVERSE   /**< The numbers of items per process are mirrored. */
};

/** All distributions of \ref t8_test_partition_shape. */
static const t8_test_partition_shape t8_test_partition_shapes[]
  = { T8_TEST_PARTITION_SHIFT, T8_TEST_PARTITION_GATHER, T8_TEST_PARTITION_SCATTER, T8_TEST_PARTITION_IDENTITY,
      T8_TEST_PARTITION_REVERSE };

/**
 * Compute global offsets from the number of items per process.
 *
 * \param[in] counts The number of items of each process.
 * \return The offsets, one entry more than \a counts.
 */
static std::vector<t8_gloidx_t>
t8_test_offsets_from_counts (const std::vector<t8_gloidx_t> &counts)
{
  std::vector<t8_gloidx_t> offsets (counts.size () + 1, 0);
  std::partial_sum (counts.begin (), counts.end (), offsets.begin () + 1);
  return offsets;
}

/**
 * Compute the old and new global offsets of a test distribution. Every process computes the same offsets
 * without communication. Process p holds 10 + p items in the initial distribution.
 *
 * \param[in] shape The distribution to compute.
 * \param[in] num_ranks The number of processes.
 * \return The old and the new offsets.
 */
static std::pair<std::vector<t8_gloidx_t>, std::vector<t8_gloidx_t>>
t8_test_partition_offsets (const t8_test_partition_shape shape, const int num_ranks)
{
  std::vector<t8_gloidx_t> counts (num_ranks);
  std::vector<t8_gloidx_t> reversed_counts (num_ranks);
  for (int iproc = 0; iproc < num_ranks; ++iproc) {
    counts[iproc] = 10 + iproc;
    reversed_counts[iproc] = 10 + (num_ranks - 1 - iproc);
  }
  const std::vector<t8_gloidx_t> initial = t8_test_offsets_from_counts (counts);
  const t8_gloidx_t num_global = initial.back ();
  /* All items on process 0. */
  std::vector<t8_gloidx_t> gathered (num_ranks + 1, num_global);
  gathered[0] = 0;

  switch (shape) {
  case T8_TEST_PARTITION_SHIFT: {
    std::vector<t8_gloidx_t> shifted = initial;
    for (int iproc = 1; iproc < num_ranks; ++iproc) {
      shifted[iproc] = std::min (initial[iproc] + 3, num_global);
    }
    return { initial, shifted };
  }
  case T8_TEST_PARTITION_GATHER:
    return { initial, gathered };
  case T8_TEST_PARTITION_SCATTER:
    return { gathered, initial };
  case T8_TEST_PARTITION_IDENTITY:
    return { initial, initial };
  case T8_TEST_PARTITION_REVERSE:
    return { initial, t8_test_offsets_from_counts (reversed_counts) };
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

/** Test repartitioning with the distributions of \ref t8_test_partition_shape. */
struct data_partition_test: public testing::TestWithParam<t8_test_partition_shape>
{
 protected:
  void
  SetUp () override
  {
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    std::tie (offsets_old, offsets_new) = t8_test_partition_offsets (GetParam (), mpisize);
    pattern
      = t8_data_partition_pattern_from_offsets (offsets_old.data (), offsets_new.data (), mpirank, mpisize, comm, 0);
  }

  /**
   * Create the data of this process in a distribution: Each item holds its global index.
   *
   * \param[in] offsets The global offsets of the distribution.
   * \return The items of this process.
   */
  std::vector<double>
  global_index_data (const std::vector<t8_gloidx_t> &offsets) const
  {
    std::vector<double> data (offsets[mpirank + 1] - offsets[mpirank]);
    std::iota (data.begin (), data.end (), static_cast<double> (offsets[mpirank]));
    return data;
  }

  int mpirank;                          /**< The rank of this process. */
  int mpisize;                          /**< The number of processes. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD; /**< The communicator of the test. */
  std::vector<t8_gloidx_t> offsets_old; /**< The global offsets before the repartitioning. */
  std::vector<t8_gloidx_t> offsets_new; /**< The global offsets after the repartitioning. */
  t8_data_partition_pattern pattern;    /**< The partition pattern of this process. */
};

/**
 * The chunks of the pattern cover the old and new local items exactly once, without self and in rank order.
 */
TEST_P (data_partition_test, pattern)
{
  EXPECT_EQ (pattern.num_old, static_cast<size_t> (offsets_old[mpirank + 1] - offsets_old[mpirank]));
  EXPECT_EQ (pattern.num_new, static_cast<size_t> (offsets_new[mpirank + 1] - offsets_new[mpirank]));

  size_t num_sent = pattern.self_count;
  size_t next_first = 0;
  int previous_rank = -1;
  for (const auto &chunk : pattern.sends) {
    EXPECT_NE (chunk.rank, mpirank);
    EXPECT_GT (chunk.rank, previous_rank);
    EXPECT_GT (chunk.count, static_cast<size_t> (0));
    /* The global indices of the chunk belong to chunk.rank in the new distribution. */
    const t8_gloidx_t global_first = offsets_old[mpirank] + static_cast<t8_gloidx_t> (chunk.first);
    EXPECT_LE (offsets_new[chunk.rank], global_first);
    EXPECT_LE (global_first + static_cast<t8_gloidx_t> (chunk.count), offsets_new[chunk.rank + 1]);
    EXPECT_GE (chunk.first, next_first);
    next_first = chunk.first + chunk.count;
    previous_rank = chunk.rank;
    num_sent += chunk.count;
  }
  EXPECT_EQ (num_sent, pattern.num_old);

  size_t num_received = pattern.self_count;
  previous_rank = -1;
  for (const auto &chunk : pattern.recvs) {
    EXPECT_NE (chunk.rank, mpirank);
    EXPECT_GT (chunk.rank, previous_rank);
    EXPECT_GT (chunk.count, static_cast<size_t> (0));
    /* The global indices of the chunk belong to chunk.rank in the old distribution. */
    const t8_gloidx_t global_first = offsets_new[mpirank] + static_cast<t8_gloidx_t> (chunk.first);
    EXPECT_LE (offsets_old[chunk.rank], global_first);
    EXPECT_LE (global_first + static_cast<t8_gloidx_t> (chunk.count), offsets_old[chunk.rank + 1]);
    previous_rank = chunk.rank;
    num_received += chunk.count;
  }
  EXPECT_EQ (num_received, pattern.num_new);

  /* The staying items have the same global indices before and after. */
  if (pattern.self_count > 0) {
    EXPECT_EQ (offsets_old[mpirank] + static_cast<t8_gloidx_t> (pattern.self_first_old),
               offsets_new[mpirank] + static_cast<t8_gloidx_t> (pattern.self_first_new));
  }
}

/**
 * Repartition plain bytes, with the default message size and with messages split at odd byte counts.
 */
TEST_P (data_partition_test, bytes)
{
  const std::vector<double> data_in = global_index_data (offsets_old);
  const std::vector<double> expected = global_index_data (offsets_new);

  for (const size_t max_message_bytes : { static_cast<size_t> (std::numeric_limits<int>::max ()), size_t { 7 } }) {
    std::vector<double> data_out (pattern.num_new, -1.0);
    t8_data_partition_bytes (pattern, data_in.data (), data_out.data (), sizeof (double), max_message_bytes);
    EXPECT_EQ (data_out, expected) << "max_message_bytes " << max_message_bytes;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_data_exchange, data_partition_test, testing::ValuesIn (t8_test_partition_shapes));

/** The number of values of every \ref t8_test_cell in the handler tests. */
static constexpr int t8_test_cell_num_values = 3;

/**
 * Create a test item that encodes its global index.
 * pseudo_tree items have a size that varies with the index and contain a nested handler for some indices.
 *
 * \tparam TType The type of the item: double, enlarged_data<int>, pseudo_tree or t8_test_cell.
 * \param[in] global_index The global index of the item.
 * \return The item.
 */
template <typename TType>
static TType
t8_test_item_from_global_index (const t8_gloidx_t global_index)
{
  if constexpr (std::is_same_v<TType, double>) {
    return static_cast<double> (global_index);
  }
  else if constexpr (std::is_same_v<TType, enlarged_data<int>>) {
    return enlarged_data<int> (static_cast<int> (global_index), static_cast<int> (2 * global_index));
  }
  else if constexpr (std::is_same_v<TType, t8_test_cell>) {
    t8_test_cell cell;
    for (int ivalue = 0; ivalue < t8_test_cell_num_values; ++ivalue) {
      cell.values.push_back (global_index + 0.25 * ivalue);
    }
    return cell;
  }
  else {
    static_assert (std::is_same_v<TType, pseudo_tree>);
    pseudo_tree tree;
    tree.topo_data.assign (global_index % 4, static_cast<int> (global_index));
    if (global_index % 2 == 1) {
      std::vector<enlarged_data<double>> values;
      for (t8_gloidx_t ivalue = 0; ivalue < global_index % 3; ++ivalue) {
        values.emplace_back (global_index + 0.5, static_cast<int> (ivalue));
      }
      tree.tree_data.push_back (std::make_shared<t8_vector_handler<enlarged_data<double>>> (std::move (values)));
    }
    return tree;
  }
}

/**
 * Check that an item encodes a given global index, see \ref t8_test_item_from_global_index.
 *
 * \tparam TType The type of the item.
 * \param[in] item The item to check.
 * \param[in] global_index The expected global index.
 * \return True if \a item encodes \a global_index, false otherwise.
 */
template <typename TType>
static bool
t8_test_item_has_global_index (const TType &item, const t8_gloidx_t global_index)
{
  if constexpr (std::is_same_v<TType, enlarged_data<int>>) {
    return item.data == global_index && item.check == 2 * global_index;
  }
  else if constexpr (std::is_same_v<TType, pseudo_tree>) {
    const pseudo_tree expected = t8_test_item_from_global_index<pseudo_tree> (global_index);
    if (item.topo_data != expected.topo_data || item.tree_data.size () != expected.tree_data.size ()) {
      return false;
    }
    if (item.tree_data.empty ()) {
      return true;
    }
    const auto *values = dynamic_cast<const t8_vector_handler<enlarged_data<double>> *> (item.tree_data[0].get ());
    const auto *expected_values
      = dynamic_cast<const t8_vector_handler<enlarged_data<double>> *> (expected.tree_data[0].get ());
    if (values == nullptr || values->get_data ().size () != expected_values->get_data ().size ()) {
      return false;
    }
    for (size_t ivalue = 0; ivalue < values->get_data ().size (); ++ivalue) {
      if (values->get_data ()[ivalue].data != expected_values->get_data ()[ivalue].data
          || values->get_data ()[ivalue].check != expected_values->get_data ()[ivalue].check) {
        return false;
      }
    }
    return true;
  }
  else {
    return item == t8_test_item_from_global_index<TType> (global_index);
  }
}

/**
 * Create a vector handler for test items. Cells get a passed data handler, all other types their t8_data_handler.
 *
 * \tparam TType The type of the items.
 * \param[in] data The items.
 * \return The handler owning \a data.
 */
template <typename TType>
static auto
t8_test_make_handler (std::vector<TType> data)
{
  if constexpr (std::is_same_v<TType, t8_test_cell>) {
    return t8_vector_handler (std::move (data), t8_test_cell_handler { t8_test_cell_num_values });
  }
  else {
    return t8_vector_handler<TType> (std::move (data));
  }
}

/**
 * Repartition vector handlers: double takes the bytewise path, enlarged_data<int> (fixed size), pseudo_tree
 * (varying size, nested handlers) and t8_test_cell (passed data handler) are serialized.
 *
 * \tparam TType The type of the items.
 */
template <typename TType>
struct handler_partition_test: public testing::Test
{
};

/** The item types of \ref handler_partition_test. */
using PartitionDataTypes = ::testing::Types<double, enlarged_data<int>, pseudo_tree, t8_test_cell>;
TYPED_TEST_SUITE (handler_partition_test, PartitionDataTypes);

/**
 * Repartition a handler with every test distribution through the abstract interface and check that every item
 * arrives at the position of its global index.
 */
TYPED_TEST (handler_partition_test, all_shapes)
{
  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int mpirank, mpisize;
  int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  for (const t8_test_partition_shape shape : t8_test_partition_shapes) {
    const auto [offsets_old, offsets_new] = t8_test_partition_offsets (shape, mpisize);
    const t8_data_partition_pattern pattern
      = t8_data_partition_pattern_from_offsets (offsets_old.data (), offsets_new.data (), mpirank, mpisize, comm, 0);

    std::vector<TypeParam> data;
    for (t8_gloidx_t iglobal = offsets_old[mpirank]; iglobal < offsets_old[mpirank + 1]; ++iglobal) {
      data.push_back (t8_test_item_from_global_index<TypeParam> (iglobal));
    }
    auto handler = t8_test_make_handler (std::move (data));
    t8_abstract_vector_handler &abstract_handler = handler;
    abstract_handler.partition (pattern);

    const std::vector<TypeParam> &new_data = handler.get_data ();
    ASSERT_EQ (new_data.size (), pattern.num_new) << "shape " << shape;
    for (size_t iitem = 0; iitem < new_data.size (); ++iitem) {
      EXPECT_TRUE (
        t8_test_item_has_global_index (new_data[iitem], offsets_new[mpirank] + static_cast<t8_gloidx_t> (iitem)))
        << "shape " << shape << ", item " << iitem;
    }
  }
}
