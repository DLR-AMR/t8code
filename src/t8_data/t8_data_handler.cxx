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

#include <test/t8_data/t8_data_handler_specs.hxx>

template <>
int
t8_data_handler<enlarged_data<int>>::t8_data_size (sc_MPI_Comm comm)
{
  int size;
  const int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &size);
  SC_CHECK_MPI (mpiret);
  return size;
}

template <>
void
t8_data_handler<enlarged_data<int>>::t8_data_pack (enlarged_data<int> data, int &pos, std::vector<char> &buffer,
                                                   sc_MPI_Comm comm)
{
  int mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Pack (&data.data, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
  SC_CHECK_MPI (mpiret);
}

template <>
void
t8_data_handler<enlarged_data<int>>::t8_data_unpack (std::vector<char> &buffer, int &pos, enlarged_data<int> &data,
                                                     int &outcount, sc_MPI_Comm comm)
{
  int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.data, 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &data.check, 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
}
