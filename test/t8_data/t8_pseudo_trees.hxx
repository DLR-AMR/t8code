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

#ifndef T8_PSEUDO_TREES_HXX
#define T8_PSEUDO_TREES_HXX

#include <t8.h>
#include <vector>
#include <algorithm>
#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler.hxx>

class pseudo_tree {
 public:
  std::vector<int> topo_data;
  std::vector<t8_abstract_data_handler *> tree_data;
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
    if (item.tree_data.size () > 0) {
      for (auto ihandler : item.tree_data) {
        total_size += ihandler->buffer_size (comm) + int_size;
      }
    }
    return total_size;
  }

  void
  pack (const pseudo_tree &data, int &pos, std::vector<char> &buffer, sc_MPI_Comm comm)
  {
    const int data_size = data.topo_data.size ();
    /* Pack number of topological data */
    int mpiret = sc_MPI_Pack (&data_size, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);
    for (int topo_item : data.topo_data) {
      /* Pack each topological data*/
      mpiret = sc_MPI_Pack (&topo_item, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
      SC_CHECK_MPI (mpiret);
    }
    /* Pack number of tree-specific data*/
    const int tree_data_size = data.tree_data.size ();
    mpiret = sc_MPI_Pack (&tree_data_size, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (auto handler : data.tree_data) {
      const int type = handler->type ();
      /* Pack type of tree data */
      mpiret = sc_MPI_Pack (&type, 1, sc_MPI_INT, buffer.data (), buffer.size (), &pos, comm);
      SC_CHECK_MPI (mpiret);
      /* Pack each data. */
      handler->pack_vector_prefix (buffer, pos, comm);
    }
  }

  void
  unpack (const std::vector<char> &buffer, int &pos, pseudo_tree &data, sc_MPI_Comm comm)
  {
    /* Unpack number of topological data */
    int topo_data_size = 0;
    int mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &topo_data_size, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    data.topo_data.resize (topo_data_size);
    for (int &topo_item : data.topo_data) {
      /* Unpack each topological item */
      mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &topo_item, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);
    }
    /* Unpack number of tree-specific data */
    int num_handler = 0;
    mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &num_handler, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    data.tree_data.resize (num_handler);

    for (auto &ihandler : data.tree_data) {
      int type;
      mpiret = sc_MPI_Unpack (buffer.data (), buffer.size (), &pos, &type, 1, sc_MPI_INT, comm);
      int outcount = 0;
      if (type == 0) {
        ihandler = new t8_data_handler<enlarged_data<int>> ();
      }
      else {
        ihandler = new t8_data_handler<enlarged_data<double>> ();
      }
      /* Unpack tree data*/
      ihandler->unpack_vector_prefix (buffer, pos, outcount, comm);
    }
  }

  int
  type ()
  {
    return -1;
  }
};

#endif /* T8_PSEUDO_TREES_HXX */
