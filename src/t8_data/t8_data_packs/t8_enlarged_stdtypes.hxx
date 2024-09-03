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

#ifndef T8_ENLARGED_STDTYPES
#define T8_ENLARGED_STDTYPES

#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler_base.hxx>

enum pseudo_types { T8_ENLARGED_INT = 0, T8_ENLARGED_DOUBLE = 1 };

template <>
class t8_single_data_handler<enlarged_data<int>> {
 public:
  int
  size (const enlarged_data<int> &item, sc_MPI_Comm comm)
  {
    int size;
    const int mpiret = sc_MPI_Pack_size (2, sc_MPI_INT, comm, &size);
    SC_CHECK_MPI (mpiret);
    return size;
  }

  void
  pack (const enlarged_data<int> &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&(data.data), 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&(data.check), 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, enlarged_data<int> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &(data.data), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &(data.check), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }

  int
  type ()
  {
    return T8_ENLARGED_INT;
  }

  void *
  new_user_handler (const int type)
  {
    return new t8_data_handler<enlarged_data<int>> ();
  }
};

template <>
class t8_single_data_handler<enlarged_data<double>> {
 public:
  int
  size (const enlarged_data<double> &item, sc_MPI_Comm comm)
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
  pack (const enlarged_data<double> &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Pack (&(data.data), 1, sc_MPI_DOUBLE, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Pack (&(data.check), 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, enlarged_data<double> &data, sc_MPI_Comm comm)
  {
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &(data.data), 1, sc_MPI_DOUBLE, comm);
    SC_CHECK_MPI (mpiret);

    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &(data.check), 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
  }

  int
  type ()
  {
    return T8_ENLARGED_DOUBLE;
  }

  void *
  new_user_handler (const int type)
  {
    return new t8_data_handler<enlarged_data<double>> ();
  }
};

template <>
class t8_single_data_handler<pseudo_tree> {
 public:
  int
  size (const pseudo_tree &item, sc_MPI_Comm comm)
  {
    int int_size = 0;
    const int topo_data_size = item.topo_data.size () + 1;

    const int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &int_size);
    SC_CHECK_MPI (mpiret);
    int total_size = topo_data_size * int_size;

    /* tree_data_size */
    total_size += int_size;

    for (t8_abstract_data_handler *handler : item.tree_data) {
      /* type_size */
      total_size += int_size;
      /* Extra int to send type of the data. */
      total_size += handler->buffer_size (comm);
    }
    return total_size;
  }

  void
  pack (const pseudo_tree &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    int data_size = data.topo_data.size ();
    int mpiret = sc_MPI_Pack (&data_size, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
    for (int topo_item : data.topo_data) {
      mpiret = sc_MPI_Pack (&topo_item, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
      SC_CHECK_MPI (mpiret);
    }
    int tree_data_size = data.tree_data.size ();
    mpiret = sc_MPI_Pack (&tree_data_size, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (t8_abstract_data_handler *handler : data.tree_data) {
      int type = handler->type ();
      mpiret = sc_MPI_Pack (&type, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
      SC_CHECK_MPI (mpiret);
      handler->pack_vector_prefix (buffer, pos, comm);
    }
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, pseudo_tree &data, sc_MPI_Comm comm)
  {
    int topo_data_size = 0;
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &topo_data_size, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    data.topo_data.resize (topo_data_size);
    for (int &topo_item : data.topo_data) {
      mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &topo_item, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);
    }

    int num_handler = 0;
    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &num_handler, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    data.tree_data.resize (num_handler);

    for (t8_abstract_data_handler *handler : data.tree_data) {
      int type;
      mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &type, 1, sc_MPI_INT, comm);
      handler = (t8_abstract_data_handler *) handler->new_handler (type);
      int outcount = 0;
      handler->unpack_vector_prefix (buffer, pos, outcount, comm);
    }
  }

  int
  type ()
  {
    return -1;
  }

  void *
  new_user_handler (const int type)
  {
    switch (type) {
    case T8_ENLARGED_INT:
      return new t8_data_handler<enlarged_data<int>> ();
    case T8_ENLARGED_DOUBLE:
      return new t8_data_handler<enlarged_data<double>> ();
    default:
      SC_ABORT_NOT_REACHED ();
      break;
    }
  }
};

#endif /* T8_ENLARGED_STDTYPES */
