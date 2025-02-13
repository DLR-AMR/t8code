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

#include <gtest/gtest.h>
#include <test/t8_data/t8_enlarged_stdtypes.hxx>
#include <test/t8_data/t8_data_handler_specs.hxx>
#include <test/t8_data/t8_pseudo_trees.hxx>
#include <t8_data/t8_data_handler.hxx>
#include <vector>
#include <numeric>

/**
 * \file Test to check the functionality of the t8_data_handler class. 
 * 
 */

/**
 * Templated testing class. Creates enlarged data (original data + a checking integer) and a 
 * data handler. 
 * 
 * @tparam T the type of data
 */
template <typename TType>
class data_handler_test: public testing::Test {
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
    data_handler = new t8_data_handler<enlarged_data<TType>> (creator.large_data);
  }

  void
  TearDown () override
  {
    delete data_handler;
  }

  t8_data_handler<enlarged_data<TType>> *data_handler;
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

  this->recv_data = *(this->data_handler->get_data ());

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
  int mpiret = this->data_handler->send (send_to, 0, this->comm);
  SC_CHECK_MPI (mpiret);

  /* Compute the rank we this rank receives from. */
  int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpirank - 1);

  /* Receive and unpack the data. */
  sc_MPI_Status status;
  int outcount;
  mpiret = this->data_handler->recv (recv_from, 0, this->comm, &status, outcount);
  SC_CHECK_MPI (mpiret);
  this->recv_data = *(this->data_handler->get_data ());

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
  t8_data_handler<enlarged_data<int>> int_handler (int_data);
  t8_data_handler<enlarged_data<double>> double_handler (double_data);
  std::vector<t8_abstract_data_handler *> handler;

  handler.push_back (&int_handler);
  handler.push_back (&double_handler);

#if T8_ENABLE_MPI
  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (mpirank + 1) % mpisize;
  /* Compute the rank this rank receives from. */
  int recv_from = (mpirank == 0) ? (mpisize - 1) : (mpirank - 1);
  for (t8_abstract_data_handler *ihandler : handler) {

    mpiret = ihandler->send (send_to, 0, comm);
    SC_CHECK_MPI (mpiret);
    /* Receive and unpack the data. */
    sc_MPI_Status status;
    int outcount;
    mpiret = ihandler->recv (recv_from, 0, comm, &status, outcount);
  }

  std::vector<enlarged_data<int>> recv_ints = *((t8_data_handler<enlarged_data<int>> *) (handler[0]))->get_data ();
  std::vector<enlarged_data<double>> recv_doubles
    = *((t8_data_handler<enlarged_data<double>> *) (handler[1]))->get_data ();

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
  tree.tree_data[0] = std::make_shared<t8_data_handler<enlarged_data<int>>> (std::move (int_data));

  pseudo_tree tree_copy (tree);
  EXPECT_EQ (tree.topo_data.size (), tree_copy.topo_data.size ());
  EXPECT_EQ (tree.tree_data.size (), tree_copy.tree_data.size ());

  auto handler = std::dynamic_pointer_cast<t8_data_handler<enlarged_data<int>>> (tree_copy.tree_data[0]).get ();
  ASSERT_NE (handler, nullptr);
  std::vector<enlarged_data<int>> copied_data = *(handler->get_data ());

  for (int idata = 0; idata < num_data; ++idata) {
    EXPECT_EQ (copied_data[idata].data, idata);
    EXPECT_EQ (copied_data[idata].check, 42);
  }

  pseudo_tree tree_equal = std::move (tree_copy);

  EXPECT_EQ (tree.topo_data.size (), tree_equal.topo_data.size ());
  EXPECT_EQ (tree.tree_data.size (), tree_equal.tree_data.size ());

  auto handler_equal = std::dynamic_pointer_cast<t8_data_handler<enlarged_data<int>>> (tree_equal.tree_data[0]).get ();
  ASSERT_NE (handler_equal, nullptr);
  std::vector<enlarged_data<int>> equal_data = *(handler_equal->get_data ());

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
        tree.tree_data[itree_data] = std::make_shared<t8_data_handler<enlarged_data<int>>> (std::move (int_data));
      }
      else {
        std::vector<enlarged_data<double>> double_data (num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          double_data[idata].check = mpirank;
          double_data[idata].data = static_cast<double> (idata) + fraction;
        }
        tree.tree_data[itree_data] = std::make_shared<t8_data_handler<enlarged_data<double>>> (std::move (double_data));
      }
    }
    trees[itree] = std::move (tree);
  }

  t8_data_handler<pseudo_tree> tree_handler (trees);

  const int send_to = (mpirank + 1) % mpisize;
  const int recv_from = (mpirank == 0) ? (mpisize - 1) : (mpirank - 1);

  mpiret = tree_handler.send (send_to, 0, comm);
  SC_CHECK_MPI (mpiret);
  sc_MPI_Status status;
  int outcount;
  mpiret = tree_handler.recv (recv_from, 0, comm, &status, outcount);

  std::vector<pseudo_tree> recv_trees = *(tree_handler.get_data ());

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
          = std::dynamic_pointer_cast<t8_data_handler<enlarged_data<int>>> (recv_trees[itree].tree_data[itree_data])
              .get ();
        ASSERT_NE (int_handler, nullptr);
        std::vector<enlarged_data<int>> recv_ints = *(int_handler->get_data ());

        ASSERT_EQ (static_cast<int> (recv_ints.size ()), num_data);
        for (int idata = 0; idata < num_data; ++idata) {
          EXPECT_EQ (recv_ints[idata].data, idata);
          EXPECT_EQ (recv_ints[idata].check, recv_from);
        }
      }
      else {
        auto double_handler
          = std::dynamic_pointer_cast<t8_data_handler<enlarged_data<double>>> (recv_trees[itree].tree_data[itree_data])
              .get ();
        ASSERT_NE (double_handler, nullptr);
        std::vector<enlarged_data<double>> recv_double = *(double_handler->get_data ());
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
