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
#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler.hxx>
#include <vector>

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

    data_handler = new t8_data_handler<enlarged_data<T>> ();
    creator = new data_creator<enlarged_data<T>> ();
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
 * Test to pack and unpack a single element of type T. 
 */
TYPED_TEST_P (data_handler_test, pack_unpack_single_data)
{
  /* Send buffer to be filled with the packed data. */
  std::vector<char> buffer (this->data_handler->buffer_size (1, this->comm));

  /* Create enlarged data. */
  this->creator->create (1);
  int pos = 0;

  /* Pack the data into the buffer. */
  this->data_handler->data_pack (this->creator->large_data[0], pos, buffer, this->comm);

  /* Unpack the data. */
  this->recv_data.resize (1);
  pos = 0;
  this->data_handler->data_unpack (buffer, pos, this->recv_data[0], this->comm);

  EXPECT_EQ (this->recv_data[0].data, this->creator->large_data[0].data);
  EXPECT_EQ (this->recv_data[0].check, this->creator->large_data[0].check);
}

/**
 * Test to pack and unpack a vector of elements of type T. 
 */
TYPED_TEST_P (data_handler_test, pack_unpack_vector_of_data)
{
  /* Test different sizes. */
  for (int num_data = 1; num_data < this->max_num_data; num_data++) {
    this->creator->create (num_data);

    /* Create send buffer and pack data into it. */
    std::vector<char> buffer (this->data_handler->buffer_size (num_data, this->comm));
    this->data_handler->data_pack_vector (this->creator->large_data, buffer, this->comm);

    int outcount = 0;
    this->data_handler->data_unpack_vector (buffer, this->recv_data, outcount, this->comm);
    EXPECT_EQ (outcount, num_data);
    for (int idata = 0; idata < num_data; idata++) {
      EXPECT_EQ (this->recv_data[idata].data, this->creator->large_data[idata].data);
      EXPECT_EQ (this->recv_data[idata].check, this->creator->large_data[idata].check);
    }
  }
}

/**
 * Use the send and receive routines for packed data. 
 */
TYPED_TEST_P (data_handler_test, send_recv)
{
  /* Create the data. */
  this->creator->create (this->max_num_data);

  /* Compute the rank this rank sends to. We send in a round-robin fashion */
  int send_to = (this->mpirank + 1) % this->mpisize;

  /* Pack and send the data. */
  int mpiret = this->data_handler->data_send (this->creator->large_data, send_to, 0, this->comm);
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
  mpiret = this->data_handler->data_recv (this->recv_data, recv_from, 0, this->comm, &status, outcount);
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

TYPED_TEST_P (data_handler_test, allgather)
{
  this->creator->create (this->max_num_data);

  const int mpiret = this->data_handler->allgather (this->creator->large_data, this->recv_data, this->comm);
  SC_CHECK_MPI (mpiret);

  const int data_size = this->mpisize * this->max_num_data;
  for (int idata = 0; idata < data_size; ++idata) {
    EXPECT_EQ (this->recv_data[idata].data, this->creator->large_data[idata % this->max_num_data].data);
    EXPECT_EQ (this->recv_data[idata].check, this->creator->large_data[idata % this->max_num_data].check);
  }
}

REGISTER_TYPED_TEST_SUITE_P (data_handler_test, pack_unpack_single_data, pack_unpack_vector_of_data, send_recv,
                             allgather);

using DataTypes = ::testing::Types<int, double>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_data_handler, data_handler_test, DataTypes, );
