/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

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
 * \file t8_gtest_data_handler.cxx
 *
 * Test to check the functionality of the t8_data_handler class.
 *
 */

#include <gtest/gtest.h>
#include <t8_data/t8_enlarged_stdtypes.hxx>
#include <t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_pseudo_trees.hxx>
#include <t8_data/t8_vector_handler.hxx>
#include <vector>
#include <numeric>
#include <type_traits>
#include <utility>

/**
 * Templated testing class. Creates enlarged data (original data + a checking integer) and a
 * data handler.
 *
 * \tparam TType the type of data
 */
template <typename TType>
struct data_handler_test: public testing::Test
{
 protected:
  void
  SetUp () override
  {
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    creator = data_creator<enlarged_data<TType>> ();
    creator.create (max_num_data);
    data_handler = new t8_vector_handler<enlarged_data<TType>> (creator.large_data);
  }

  void
  TearDown () override
  {
    delete data_handler;
  }

  t8_vector_handler<enlarged_data<TType>> *data_handler;
  data_creator<enlarged_data<TType>> creator;
  std::vector<enlarged_data<TType>> recv_data;
  int mpirank;
  int mpisize;
  const int max_num_data = 100;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

TYPED_TEST_SUITE_P (data_handler_test);

/**
 * Test to pack and unpack a vector of elements of type T.
 */
TYPED_TEST_P (data_handler_test, pack_unpack_vector_of_data)
{
  /* Create send buffer and pack data into it. */
  int pos = 0;
  const int num_bytes = this->data_handler->buffer_size (this->comm);
  void *buffer = malloc (num_bytes);
  this->data_handler->pack_vector_prefix (buffer, num_bytes, pos, this->comm);

  int outcount = 0;
  pos = 0;
  this->data_handler->unpack_vector_prefix (buffer, num_bytes, pos, outcount, this->comm);
  EXPECT_EQ (outcount, this->max_num_data);

  this->recv_data = this->data_handler->get_data ();

  for (int idata = 0; idata < this->max_num_data; idata++) {
    EXPECT_EQ (this->recv_data[idata].data, this->creator.large_data[idata].data);
    EXPECT_EQ (this->recv_data[idata].check, this->creator.large_data[idata].check);
  }

  free (buffer);
}

/**
 * Use the send and receive routines for packed data.
 */
TYPED_TEST_P (data_handler_test, send_recv)
{
  /* Compute the rank this rank sends to. We send in a round-robin fashion */

  /* Pack and send the data. */
#if T8_ENABLE_MPI
  int send_to = (this->mpirank + 1) % this->mpisize;
  if (send_to == this->mpirank) {
    t8_debugf ("Rank %d not sending data to itself\n", this->mpirank);
    return;
  }
  this->data_handler->send (send_to, 0, this->comm);

  /* Compute the rank we this rank receives from. */
  int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpirank - 1);

  /* Receive and unpack the data. */
  sc_MPI_Status status;
  int outcount;
  this->data_handler->recv (recv_from, 0, this->comm, &status, outcount);
  this->recv_data = this->data_handler->get_data ();

  EXPECT_EQ (outcount, this->max_num_data);
  for (int idata = 0; idata < this->max_num_data; idata++) {
    EXPECT_EQ (this->recv_data[idata].data, this->creator.large_data[idata].data);
    EXPECT_EQ (this->recv_data[idata].check, this->creator.large_data[idata].check);
  }
#endif
}

TEST (data_handler_test, multiple_handler)
{
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int mpirank;
  int mpisize;
  int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  const int num_data = 1000;

  std::vector<enlarged_data<int>> int_data (num_data);
  std::vector<enlarged_data<double>> double_data (num_data);
  const double fraction = 0.42;
  for (int idata = 0; idata < num_data; idata++) {
    int_data[idata].data = idata;
    int_data[idata].check = mpirank;
    double_data[idata].data = (double) idata + fraction;
    double_data[idata].check = mpirank;
  }
  t8_vector_handler<enlarged_data<int>> int_handler (int_data);
  t8_vector_handler<enlarged_data<double>> double_handler (double_data);

#if T8_ENABLE_MPI
  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (mpirank + 1) % mpisize;
  /* Compute the rank this rank receives from. */
  int recv_from = (mpirank - 1 + mpisize) % mpisize;
  if (send_to != mpirank) {
    sc_MPI_Status status;
    int outcount;
    int_handler.sendrecv (send_to, recv_from, 0, comm, &status, outcount);
    double_handler.sendrecv (send_to, recv_from, 0, comm, &status, outcount);
  }
  else {
    t8_debugf ("Rank %i not sending data to itself\n", mpirank);
  }

  const auto &recv_ints = int_handler.get_data ();
  const auto &recv_doubles = double_handler.get_data ();

  SC_CHECK_MPI (mpiret);
  for (int idata = 0; idata < num_data; idata++) {
    EXPECT_EQ (recv_ints[idata].check, recv_from);
    EXPECT_EQ (recv_ints[idata].data, idata);
    EXPECT_EQ (recv_doubles[idata].check, recv_from);
    EXPECT_NEAR (recv_doubles[idata].data, (double) idata + fraction, T8_PRECISION_EPS);
  }
#endif
}

TEST (data_handler_test, pseudo_tree_test)
{
  const int num_data = 100;
  pseudo_tree tree;
  tree.topo_data.resize (10);
  std::iota (tree.topo_data.begin (), tree.topo_data.end (), 0);

  std::vector<enlarged_data<int>> int_data (num_data);
  for (int idata = 0; idata < num_data; ++idata) {
    int_data[idata].check = 42;
    int_data[idata].data = idata;
  }
  tree.tree_data.resize (1);
  tree.tree_data[0] = std::make_shared<t8_vector_handler<enlarged_data<int>>> (std::move (int_data));

  pseudo_tree tree_copy (tree);
  EXPECT_EQ (tree.topo_data.size (), tree_copy.topo_data.size ());
  EXPECT_EQ (tree.tree_data.size (), tree_copy.tree_data.size ());

  auto handler = std::dynamic_pointer_cast<t8_vector_handler<enlarged_data<int>>> (tree_copy.tree_data[0]).get ();
  ASSERT_NE (handler, nullptr);
  const std::vector<enlarged_data<int>> &copied_data = handler->get_data ();

  for (int idata = 0; idata < num_data; ++idata) {
    EXPECT_EQ (copied_data[idata].data, idata);
    EXPECT_EQ (copied_data[idata].check, 42);
  }

  pseudo_tree tree_equal = std::move (tree_copy);

  EXPECT_EQ (tree.topo_data.size (), tree_equal.topo_data.size ());
  EXPECT_EQ (tree.tree_data.size (), tree_equal.tree_data.size ());

  auto handler_equal
    = std::dynamic_pointer_cast<t8_vector_handler<enlarged_data<int>>> (tree_equal.tree_data[0]).get ();
  ASSERT_NE (handler_equal, nullptr);
  const std::vector<enlarged_data<int>> &equal_data = handler_equal->get_data ();

  for (int idata = 0; idata < num_data; ++idata) {
    EXPECT_EQ (equal_data[idata].data, idata);
    EXPECT_EQ (equal_data[idata].check, 42);
  }
}

#if T8_ENABLE_MPI
TEST (data_handler_test, tree_test)
{
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int mpirank, mpisize, mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  const int num_trees = (mpirank % 4) * 10;
  const int num_data = 1000;
  const double fraction = 0.42;

  std::vector<pseudo_tree> trees (num_trees);

  for (int itree = 0; itree < num_trees; ++itree) {
    pseudo_tree tree;
    const int tree_topo_size = ((mpirank % 3) + 1) * 10;
    tree.topo_data.resize (tree_topo_size);
    std::iota (tree.topo_data.begin (), tree.topo_data.end (), 0);

    const int num_tree_data = (mpirank + itree) % 2;
    tree.tree_data.resize (num_tree_data);
    for (int itree_data = 0; itree_data < num_tree_data; ++itree_data) {
      if (itree_data == 0) {
        std::vector<enlarged_data<int>> int_data (num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          int_data[idata].check = mpirank;
          int_data[idata].data = idata;
        }
        tree.tree_data[itree_data] = std::make_shared<t8_vector_handler<enlarged_data<int>>> (std::move (int_data));
      }
      else {
        std::vector<enlarged_data<double>> double_data (num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          double_data[idata].check = mpirank;
          double_data[idata].data = static_cast<double> (idata) + fraction;
        }
        tree.tree_data[itree_data]
          = std::make_shared<t8_vector_handler<enlarged_data<double>>> (std::move (double_data));
      }
    }
    trees[itree] = std::move (tree);
  }

  t8_vector_handler<pseudo_tree> tree_handler (trees);

  const int send_to = (mpirank + 1) % mpisize;
  const int recv_from = (mpirank == 0) ? (mpisize - 1) : (mpirank - 1);
  if (send_to == mpirank) {
    t8_debugf ("Rank %d not sending data to itself\n", mpirank);
    return;
  }
  tree_handler.send (send_to, 0, comm);
  sc_MPI_Status status;
  int outcount;
  tree_handler.recv (recv_from, 0, comm, &status, outcount);

  const std::vector<pseudo_tree> &recv_trees = tree_handler.get_data ();

  const int num_recv_trees = recv_trees.size ();
  ASSERT_EQ (num_recv_trees, (recv_from % 4) * 10);

  for (int itree = 0; itree < num_recv_trees; ++itree) {
    const int num_recv_tree_topo_size = recv_trees[itree].topo_data.size ();
    ASSERT_EQ (num_recv_tree_topo_size, ((recv_from % 3) + 1) * 10);

    for (int itopo_data = 0; itopo_data < num_recv_tree_topo_size; ++itopo_data) {
      EXPECT_EQ (recv_trees[itree].topo_data[itopo_data], itopo_data);
    }

    const int num_recv_tree_data = recv_trees[itree].tree_data.size ();
    ASSERT_EQ (num_recv_tree_data, (recv_from + itree) % 2);
    for (int itree_data = 0; itree_data < num_recv_tree_data; ++itree_data) {
      if (itree_data == 0) {
        auto int_handler
          = std::dynamic_pointer_cast<t8_vector_handler<enlarged_data<int>>> (recv_trees[itree].tree_data[itree_data])
              .get ();
        ASSERT_NE (int_handler, nullptr);
        const std::vector<enlarged_data<int>> &recv_ints = int_handler->get_data ();

        ASSERT_EQ (static_cast<int> (recv_ints.size ()), num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          EXPECT_EQ (recv_ints[idata].data, idata);
          EXPECT_EQ (recv_ints[idata].check, recv_from);
        }
      }
      else {
        auto double_handler = std::dynamic_pointer_cast<t8_vector_handler<enlarged_data<double>>> (
                                recv_trees[itree].tree_data[itree_data])
                                .get ();
        ASSERT_NE (double_handler, nullptr);
        const std::vector<enlarged_data<double>> &recv_double = double_handler->get_data ();
        ASSERT_EQ (static_cast<int> (recv_double.size ()), num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          EXPECT_EQ (recv_double[idata].data, static_cast<double> (idata) + fraction);
          EXPECT_EQ (recv_double[idata].check, recv_from);
        }
      }
    }
  }
}
#endif

REGISTER_TYPED_TEST_SUITE_P (data_handler_test, pack_unpack_vector_of_data, send_recv);

using DataTypes = ::testing::Types<int, double>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_data_handler, data_handler_test, DataTypes, );

/**
 * A trivially copyable struct with padding. It has no t8_data_handler specialization and is therefore
 * handled by the generic handler.
 */
struct t8_test_pod
{
  int index;    /**< Index of the item. */
  char tag;     /**< A character derived from the rank. */
  double value; /**< A floating point value derived from the index. */

  /** Compare two items field by field. */
  bool
  operator== (const t8_test_pod &other) const
    = default;
};

/* Trivially copyable types without explicit specialization are communicated as raw bytes. */
static_assert (T8BytewiseData<int>);
static_assert (T8BytewiseData<double>);
static_assert (T8BytewiseData<t8_test_pod>);
/* enlarged_data is trivially copyable, but its explicit specialization packs field-wise. */
static_assert (std::is_trivially_copyable_v<enlarged_data<int>>);
static_assert (!T8BytewiseData<enlarged_data<int>>);
static_assert (!T8BytewiseData<enlarged_data<double>>);
static_assert (!T8BytewiseData<pseudo_tree>);

/**
 * Create a test item of a generic type.
 *
 * \tparam TType The type of the item.
 * \param[in] index The index of the item.
 * \param[in] rank The rank that creates the item.
 * \return An item that depends on \a index and \a rank.
 */
template <typename TType>
static TType
t8_test_generic_item (const int index, const int rank)
{
  if constexpr (std::is_same_v<TType, t8_test_pod>) {
    return t8_test_pod { index, static_cast<char> ('a' + rank % 26), index + 0.5 };
  }
  else {
    return static_cast<TType> (index + 1000 * rank);
  }
}

/**
 * Create the test data of a rank. The number of items depends on the rank.
 *
 * \tparam TType The type of the items.
 * \param[in] rank The rank that creates the data.
 * \return The test data of \a rank.
 */
template <typename TType>
static std::vector<TType>
t8_test_generic_data (const int rank)
{
  std::vector<TType> data (10 + rank);
  for (size_t iitem = 0; iitem < data.size (); ++iitem) {
    data[iitem] = t8_test_generic_item<TType> (static_cast<int> (iitem), rank);
  }
  return data;
}

/**
 * Test the generic handler of trivially copyable types without explicit specialization.
 *
 * \tparam TType The type of the data.
 */
template <typename TType>
struct generic_data_handler_test: public testing::Test
{
 protected:
  void
  SetUp () override
  {
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
  }

  int mpirank;                          /**< The rank of this process. */
  int mpisize;                          /**< The number of processes. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD; /**< The communicator of the test. */
};

using GenericDataTypes = ::testing::Types<int, double, t8_test_pod>;
TYPED_TEST_SUITE (generic_data_handler_test, GenericDataTypes);

/**
 * Pack a vector of generic items and unpack it into another handler.
 */
TYPED_TEST (generic_data_handler_test, pack_unpack)
{
  const std::vector<TypeParam> data = t8_test_generic_data<TypeParam> (this->mpirank);
  t8_vector_handler<TypeParam> handler (data);

  const int num_bytes = handler.buffer_size (this->comm);
  std::vector<char> buffer (num_bytes);
  int position = 0;
  handler.pack_vector_prefix (buffer.data (), num_bytes, position, this->comm);
  EXPECT_LE (position, num_bytes);

  t8_vector_handler<TypeParam> unpacked;
  int outcount = 0;
  position = 0;
  unpacked.unpack_vector_prefix (buffer.data (), num_bytes, position, outcount, this->comm);
  EXPECT_EQ (outcount, static_cast<int> (data.size ()));
  EXPECT_EQ (unpacked.get_data (), data);
}

/**
 * Exchange vectors of generic items of different sizes in a ring.
 */
TYPED_TEST (generic_data_handler_test, sendrecv)
{
#if T8_ENABLE_MPI
  const int send_to = (this->mpirank + 1) % this->mpisize;
  const int recv_from = (this->mpirank - 1 + this->mpisize) % this->mpisize;
  t8_vector_handler<TypeParam> handler (t8_test_generic_data<TypeParam> (this->mpirank));

  sc_MPI_Status status;
  int outcount = 0;
  handler.sendrecv (send_to, recv_from, 0, this->comm, &status, outcount);

  const std::vector<TypeParam> expected = t8_test_generic_data<TypeParam> (recv_from);
  EXPECT_EQ (outcount, static_cast<int> (expected.size ()));
  EXPECT_EQ (handler.get_data (), expected);
#endif
}

/**
 * Check the type ids of generic and explicitly specialized handlers.
 */
TEST (data_handler_test, type_ids)
{
  EXPECT_TRUE (t8_vector_handler<int> ().type () == t8_data_handler_type_unregistered);
  EXPECT_TRUE (t8_vector_handler<t8_test_pod> ().type () == t8_data_handler_type_unregistered);
  EXPECT_TRUE (t8_vector_handler<enlarged_data<int>> ().type () == t8_data_handler_type (T8_ENLARGED_INT));
  EXPECT_TRUE (t8_vector_handler<enlarged_data<double>> ().type () == t8_data_handler_type (T8_ENLARGED_DOUBLE));
}

/**
 * Data moved into and out of a handler is not copied.
 */
TEST (data_handler_test, move_set_take)
{
  std::vector<double> values (100, 1.5);
  const double *storage = values.data ();

  t8_vector_handler<double> handler (std::move (values));
  EXPECT_EQ (handler.get_data ().data (), storage);

  std::vector<double> taken = handler.take_data ();
  EXPECT_EQ (taken.data (), storage);
  EXPECT_TRUE (handler.get_data ().empty ());

  handler.set_data (std::move (taken));
  EXPECT_EQ (handler.get_data ().data (), storage);
  EXPECT_EQ (handler.get_data ().size (), static_cast<size_t> (100));
}

/* A passed data handler is accepted if it has the required members; it serializes unless it opts in to raw bytes. */
static_assert (T8DataHandler<t8_test_cell_handler, t8_test_cell>);
static_assert (!T8BytewiseDataHandler<t8_test_cell_handler, t8_test_cell>);
/* Without specialization and without passed handler, a type that is not trivially copyable has no data handler. */
static_assert (!T8DataHandler<t8_data_handler<t8_test_cell>, t8_test_cell>);

/**
 * Create the cells of a rank. The number of cells depends on the rank.
 *
 * \param[in] rank The rank that creates the cells.
 * \param[in] num_values The number of values per cell.
 * \return The cells of \a rank.
 */
static std::vector<t8_test_cell>
t8_test_cells (const int rank, const int num_values)
{
  std::vector<t8_test_cell> cells (5 + rank);
  for (size_t icell = 0; icell < cells.size (); ++icell) {
    for (int ivalue = 0; ivalue < num_values; ++ivalue) {
      cells[icell].values.push_back (1000.0 * rank + 10.0 * icell + ivalue);
    }
  }
  return cells;
}

/**
 * The data handler is deduced from the constructor arguments and its state is used for packing.
 */
TEST (data_handler_test, passed_data_handler_pack_unpack)
{
  const int num_values = 3;
  const std::vector<t8_test_cell> cells = t8_test_cells (2, num_values);
  t8_vector_handler handler (cells, t8_test_cell_handler { num_values });
  static_assert (std::is_same_v<decltype (handler), t8_vector_handler<t8_test_cell, t8_test_cell_handler>>);
  EXPECT_TRUE (handler.type () == t8_data_handler_type_unregistered);

  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  const int num_bytes = handler.buffer_size (comm);
  std::vector<char> buffer (num_bytes);
  int position = 0;
  handler.pack_vector_prefix (buffer.data (), num_bytes, position, comm);
  EXPECT_LE (position, num_bytes);

  t8_vector_handler<t8_test_cell, t8_test_cell_handler> unpacked ({}, t8_test_cell_handler { num_values });
  int outcount = 0;
  position = 0;
  unpacked.unpack_vector_prefix (buffer.data (), num_bytes, position, outcount, comm);
  EXPECT_EQ (outcount, static_cast<int> (cells.size ()));
  EXPECT_EQ (unpacked.get_data (), cells);
}

/**
 * Exchange cells with a passed data handler in a ring, with a different number of cells per rank.
 */
TEST (data_handler_test, passed_data_handler_sendrecv)
{
#if T8_ENABLE_MPI
  const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int mpirank, mpisize;
  int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  const int send_to = (mpirank + 1) % mpisize;
  const int recv_from = (mpirank - 1 + mpisize) % mpisize;

  const int num_values = 4;
  t8_vector_handler handler (t8_test_cells (mpirank, num_values), t8_test_cell_handler { num_values });
  sc_MPI_Status status;
  int outcount = 0;
  handler.sendrecv (send_to, recv_from, 0, comm, &status, outcount);

  const std::vector<t8_test_cell> expected = t8_test_cells (recv_from, num_values);
  EXPECT_EQ (outcount, static_cast<int> (expected.size ()));
  EXPECT_EQ (handler.get_data (), expected);
#endif
}
