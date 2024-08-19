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
  }

  void
  TearDown () override
  {
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

TYPED_TEST_P (data_handler_test, single_data)
{
  std::vector<char> buffer;

  this->creator->create (1);
  int pos = 0;
  this->data_handler->t8_data_pack (this->creator->large_data[0], pos, buffer, this->comm);

  int mpiret
    = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, (this->mpirank + 1) % this->mpisize, 0, this->comm);
  SC_CHECK_MPI (mpiret);

  int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpisize - 1);

  sc_MPI_Status status;

  mpiret = sc_MPI_Probe (recv_from, 0, this->comm, &status);
  SC_CHECK_MPI (mpiret);

  int size;
  mpiret = sc_MPI_Get_count (&status, sc_MPI_PACKED, &size);
  SC_CHECK_MPI (mpiret);
  std::vector<char> packed (size);

  pos = 0;
  mpiret = sc_MPI_Recv (packed.data (), packed.size (), sc_MPI_PACKED, recv_from, pos, this->comm, &status);
  SC_CHECK_MPI (mpiret);

  this->recv_data.resize (1);
  int outcount = 0;
  pos = 0;
  this->data_handler->t8_data_unpack (packed, pos, this->recv_data[0], this->comm);

  EXPECT_EQ (this->recv_data[0].data, this->creator->large_data[0].data);
  EXPECT_EQ (this->recv_data[0].check, this->creator->large_data[0].check);
}

TYPED_TEST_P (data_handler_test, vector_of_data)
{
  for (int num_data = 1; num_data < this->max_num_data; num_data++) {
    this->creator->create (num_data);

    std::vector<char> buffer;
    this->data_handler->t8_data_pack_vector (this->creator->large_data, num_data, buffer, this->comm);

    int mpiret
      = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, (this->mpirank + 1) % this->mpisize, 0, this->comm);
    SC_CHECK_MPI (mpiret);

    sc_MPI_Status status;

    int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpisize - 1);

    mpiret = sc_MPI_Probe (recv_from, 0, this->comm, &status);
    SC_CHECK_MPI (mpiret);

    int size;
    mpiret = sc_MPI_Get_count (&status, sc_MPI_PACKED, &size);
    SC_CHECK_MPI (mpiret);
    std::vector<char> packed (size);

    sc_MPI_Recv (packed.data (), packed.size (), sc_MPI_PACKED, recv_from, 0, this->comm, &status);
    SC_CHECK_MPI (mpiret);

    int outcount = 0;
    this->data_handler->t8_data_unpack_vector (packed, this->recv_data, outcount, this->comm);
    EXPECT_EQ (outcount, num_data);
    for (int idata = 0; idata < num_data; idata++) {
      EXPECT_EQ (this->recv_data[idata].data, this->creator->large_data[idata].data);
      EXPECT_EQ (this->recv_data[idata].check, this->creator->large_data[idata].check);
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P (data_handler_test, single_data, vector_of_data);

using DataTypes = ::testing::Types<int>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_data_handler, data_handler_test, DataTypes, );
