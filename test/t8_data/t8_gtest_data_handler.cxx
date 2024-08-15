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
#include <t8_data/t8_data_handler.hxx>
#include <vector>

template <typename T>
struct enlarged_data
{
  T data;     // original data
  int check;  // additional data to check against
};

template <typename T>
class data_creator {
  data_creator (const int num_data);

  std::vector<enlarged_data<T>> data;
};

template <>
class data_creator<int> {
  data_creator (const int num_data)
  {
    data.resize (num_data);

    for (int idata = 0; idata < num_data; idata++) {
      data[idata].data = 42;
      data[idata].check = idata;
    }
  }
  std::vector<enlarged_data<int>> data;
};

template <>
class data_creator<double> {
  data_creator (const int num_data)
  {
    data.resize (num_data);

    for (int idata = 0; idata < num_data; idata++) {
      data[idata].data = 42.42;
      data[idata].check = idata;
    }
  }
  std::vector<enlarged_data<int>> data;
};

template <>
class t8_data_handler<enlarged_data<int>> {
  int
  t8_data_size ()
  {
    int size;
    const int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &size);
    SC_CHECK_MPI (mpiret);
    return size;
  }

  void
  t8_data_pack (const enlarged_data<int> &data, const int pos, std::vector<char> &buffer, comm)
  {
    int size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &size);
    SC_CHECK_MPI (mpiret);
    int current_pos = pos;

    mpiret = sc_MPI_PACK (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &current_pos, comm);
    SC_CHECK_MPI (mpiret);

    current_pos += size;
    sc_MPI_PACK (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &current_pos, comm);
    SC_CHECK_MPI (mpiret);
  }
}

template <typename T>
class data_handler_test: public testing::Test {
 protected:
  void
  SetUp () override
  {
    data_handler = new t8_data_handler<T> ();

    int mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
  }

  void
  TearDown () override
  {
  }

  t8_data_handler<T> *data_handler;
  data_creator<T> creator;
  int mpirank;
  int mpisize;
  std::vector<enlarged_data<T>> recv_data;
  const int max_num_data = 100;
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
};

TYPED_TEST_SUITE_P (data_handler_test);

TYPED_TEST_P (data_handler_test, single_data)
{
  this->data_creator (1);

  std::vector<char> buffer;

  this->creator (1);
  this->data_handler->t8_data_pack (this->creator->data, 1, buffer, this->comm);

    int mpiret = sc_MPI_Send(buffer.data()
}

TYPED_TEST_P (data_handler_test, vector_of_data)
{
  for (int num_data = 1; num_data < this->max_num_data; num_data++) {
    this->creator (num_data);

    std::vector<char> buffer;
    this->data_handler->t8_data_pack_vector (this->creator.data, num_data, buffer, this->comm);

    int mpiret
      = sc_MPI_Send (buffer.data (), buffer.size (), sc_MPI_PACKED, (this->mpirank + 1) % this->mpisize, 0, this->comm);
    SC_CHECK_MPI (mpiret);

    sc_MPI_Status status;

    int recv_from = (this->mpirank == 0) ? (this->mpisize - 1) : (this->mpisize - 1);

    sc_MPI_Probe (recv_from, 0, this->comm, &status);
    SC_CHECK_MPI (mpiret);

    int size;
    sc_MPI_Get_count (&status, sc_MPI_PACKED, &size);
    SC_CHECK_MPI (mpiret);
    std::vector<char> packed (size);

    sc_MPI_Recv (packed.data (), packed.size (), sc_MPI_PACKED, recv_from, 0, this->comm, &status);
    SC_CHECK_MPI (mpiret);

    this->recv_data.resize (num_data);
    this->data_handler->t8_data_unpack (packed.data (), num_data, num_data, this->recv_data, this->comm);
    for (int idata = 0; idata < num_data; idata++) {
      EXPECT_EQ (this->recv_data[idata].data, this->creator.data[idata].data);
      EXPECT_EQ (this->recv_data[idata].check, this->creator.data[idata].check);
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P (data_handler_test, single_data, vector_of_data);

using DataTypes = ::testing::Types<int, double>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_data_handler, data_handler_test, DataTypes, );
