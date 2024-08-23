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

#ifndef T8_ENLARGED_INT
#define T8_ENLARGED_INT

#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler_base.hxx>

template <>
class t8_single_data_handler<enlarged_data<int>> {
 public:
  int
  size (sc_MPI_Comm comm)
  {
    int size;
    const int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &size);
    SC_CHECK_MPI (mpiret);
    return size;
  }

  void
  pack (const enlarged_data<int> data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, enlarged_data<int> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.data, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.check, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }
};

template <>
class t8_single_data_handler<enlarged_data<double>> {
 public:
  int
  size (sc_MPI_Comm comm)
  {
    int int_size;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &int_size);
    SC_CHECK_MPI (mpiret);
    int double_size;
    mpiret = sc_MPI_Pack_size (1, sc_MPI_DOUBLE, comm, &double_size);
    SC_CHECK_MPI (mpiret);
    return int_size + double_size;
  }

  void
  pack (const enlarged_data<double> data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_DOUBLE, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, enlarged_data<double> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.data, 1, sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.check, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }
};

#endif /* T8_ENLARGED_INT */
