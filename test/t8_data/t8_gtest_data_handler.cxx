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
#include <t8_data/t8_data_packs/t8_packed_types.hxx>
#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler.hxx>
#include <vector>
#include <numeric>

/**
 * Templated testing class. Creates enlarged data (original data + a checking integer) and a 
 * data handler. 
 * 
 * @tparam T the type of data
 */
template <typename T>
class data_handler_test: public testing::Test {
 protected:
  void
  SetUp () override
  {
    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    creator = new data_creator<enlarged_data<T>> ();
    creator->create (max_num_data);
    data_handler = new t8_data_handler<enlarged_data<T>> (creator->large_data);
  }

  void
  TearDown () override
  {
    delete (data_handler);
    delete (creator);
  }

  t8_data_handler<enlarged_data<T>> *data_handler;
  data_creator<enlarged_data<T>> *creator;
  std::vector<enlarged_data<T>> recv_data;
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
  std::vector<char> buffer (this->data_handler->buffer_size (this->comm));
  this->data_handler->pack_vector_prefix (buffer, pos, this->comm);

  int outcount = 0;
  pos = 0;
  this->data_handler->unpack_vector_prefix (buffer, pos, outcount, this->comm);
  EXPECT_EQ (outcount, this->max_num_data);

  this->recv_data = this->data_handler->get_data ();

  for (int idata = 0; idata < this->max_num_data; idata++) {
    EXPECT_EQ (this->recv_data[idata].data, this->creator->large_data[idata].data);
    EXPECT_EQ (this->recv_data[idata].check, this->creator->large_data[idata].check);
  }
}

/**
 * Use the send and receive routines for packed data. 
 */
TYPED_TEST_P (data_handler_test, send_recv)
{
  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (this->mpirank + 1) % this->mpisize;

  /* Pack and send the data. */
  int mpiret = this->data_handler->send (send_to, 0, this->comm);
#if T8_ENABLE_MPI
  SC_CHECK_MPI (mpiret);
#else
  EXPECT_EQ (mpiret, sc_MPI_ERR_OTHER);
#endif

  /* Compute the rank we this rank receives from. */
  int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpirank - 1);

  /* Receive and unpack the data. */
  sc_MPI_Status status;
  int outcount;
  mpiret = this->data_handler->recv (recv_from, 0, this->comm, &status, outcount);

  this->recv_data = this->data_handler->get_data ();
#if T8_ENABLE_MPI
  SC_CHECK_MPI (mpiret);
  EXPECT_EQ (outcount, this->max_num_data);
  for (int idata = 0; idata < this->max_num_data; idata++) {
    EXPECT_EQ (this->recv_data[idata].data, this->creator->large_data[idata].data);
    EXPECT_EQ (this->recv_data[idata].check, this->creator->large_data[idata].check);
  }
#else
  EXPECT_EQ (mpiret, sc_MPI_ERR_OTHER);
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

  t8_data_handler<enlarged_data<int>> *int_handler = new t8_data_handler<enlarged_data<int>> (int_data);
  t8_data_handler<enlarged_data<double>> *double_handler = new t8_data_handler<enlarged_data<double>> (double_data);

  std::vector<t8_abstract_data_handler *> handler = { int_handler, double_handler };

  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (mpirank + 1) % mpisize;
  /* Compute the rank this rank receives from. */
  int recv_from = (mpirank == 0) ? (mpisize - 1) : (mpirank - 1);
  t8_debugf ("[D] computed recv from %i\n", recv_from);

  for (t8_abstract_data_handler *ihandler : handler) {

    mpiret = ihandler->send (send_to, 0, comm);
#if T8_ENABLE_MPI
    SC_CHECK_MPI (mpiret);
    /* Receive and unpack the data. */
    sc_MPI_Status status;
    int outcount;
    mpiret = ihandler->recv (recv_from, 0, comm, &status, outcount);
#else
    EXPECT_EQ (mpiret, sc_MPI_ERR_OTHER);
#endif
  }

  std::vector<enlarged_data<int>> recv_ints = int_handler->get_data ();
  std::vector<enlarged_data<double>> recv_doubles = double_handler->get_data ();

#if T8_ENABLE_MPI
  SC_CHECK_MPI (mpiret);
  for (int idata = 0; idata < num_data; idata++) {
    EXPECT_EQ (recv_ints[idata].check, recv_from);
    EXPECT_EQ (recv_ints[idata].data, idata);
    EXPECT_EQ (recv_doubles[idata].check, recv_from);
    EXPECT_NEAR (recv_doubles[idata].data, (double) idata + fraction, T8_PRECISION_EPS);
  }
#else
  EXPECT_EQ (mpiret, sc_MPI_ERR_OTHER);
#endif

  delete (int_handler);
  delete (double_handler);
  /* Pack and send the data. */
}

TEST (data_handler_test, tree_test)
{
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
  int mpirank;
  int mpisize;
  int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  const int num_trees = mpirank % 4 * 10;
  const int num_data = 100;

  std::vector<pseudo_tree> trees (100);
  const int num_topo_data = ((mpirank % 2) + 1) * 100;
  std::vector<enlarged_data<int>> int_data (num_data);
  std::vector<enlarged_data<double>> double_data (num_data);
  const double fraction = 0.42;
  for (int idata = 0; idata < num_data; idata++) {
    int_data[idata].data = idata;
    int_data[idata].check = mpirank;
    double_data[idata].data = (double) idata + fraction;
    double_data[idata].check = mpirank;
  }

  for (int itree = 0; itree < num_trees; itree++) {
    std::vector<int> topo_data (num_topo_data);
    std::iota (topo_data.begin (), topo_data.end (), 0);
    const int num_tree_data = (mpirank + itree) % 3;
    std::vector<t8_abstract_data_handler *> tree_data (num_data);

    for (int idata = 0; idata < num_tree_data - 1; idata++) {
      if (idata == 0) {
        tree_data[0] = new t8_data_handler<enlarged_data<int>> (int_data);
      }
      if (idata == 1) {
        tree_data[1] = new t8_data_handler<enlarged_data<double>> (double_data);
      }
    }
    trees[itree].topo_data = topo_data;
    trees[itree].tree_data = tree_data;
  }

  t8_data_handler<pseudo_tree> *tree_handler = new t8_data_handler<pseudo_tree> (trees);

  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (mpirank + 1) % mpisize;
  int recv_from = (mpirank == 0) ? (mpisize - 1) : (mpirank - 1);

  mpiret = tree_handler->send (send_to, 0, comm);

#if T8_ENABLE_MPI
  SC_CHECK_MPI (mpiret);
#else
  EXPECT_EQ (mpiret, sc_MPI_ERR_OTHER);
#endif
  /* Receive and unpack the data. */
  sc_MPI_Status status;
  int outcount;
  tree_handler->recv (recv_from, 0, comm, &status, outcount);
}

REGISTER_TYPED_TEST_SUITE_P (data_handler_test, pack_unpack_vector_of_data, send_recv);

using DataTypes = ::testing::Types<int, double>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_data_handler, data_handler_test, DataTypes, );
