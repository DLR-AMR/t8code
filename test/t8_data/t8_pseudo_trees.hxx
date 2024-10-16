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
#include <test/t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_data_handler.hxx>
#include <memory>

/**
 *  /file This file provides a pseudo tree structure for testing purposes.
 * 
 */

/**
 * \class pseudo_tree
 * Represents a pseudo tree structure containing topological and tree data.
 * 
 * \note This class is used for testing purposes only.
 * 
 * The pseudo_tree class encapsulates two main data members:
 * - `topo_data`: A vector of integers representing topological data.
 * - `tree_data`: A vector of pointers to t8_abstract_data_handler objects, representing tree data.
 */
class pseudo_tree {
 public:
  std::vector<int> topo_data;
  std::vector<std::shared_ptr<t8_abstract_data_handler>> tree_data;
};

template <>
class t8_single_data_handler<pseudo_tree> {
 public:
  int
  size (const pseudo_tree &item, sc_MPI_Comm comm)
  {
    int int_size = 0;
    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &int_size);
    SC_CHECK_MPI (mpiret);

    // Calculate the size for topo_data
    const int topo_data_size = item.topo_data.size ();
    int total_size = (topo_data_size + 1) * int_size;

    // Calculate the size for tree_data
    total_size += int_size;  // for tree_data_size
    for (auto ihandler : item.tree_data) {
      total_size += ihandler->buffer_size (comm) + int_size;
    }

    return total_size;
  }

  void
  pack (const pseudo_tree &data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm)
  {
    const int data_size = data.topo_data.size ();
    /* Pack number of topological data */
    int mpiret = sc_MPI_Pack (&data_size, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack topological data in one call */
    mpiret = sc_MPI_Pack (data.topo_data.data (), data_size, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack number of tree-specific data */
    const int tree_data_size = data.tree_data.size ();
    mpiret = sc_MPI_Pack (&tree_data_size, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (const auto &handler : data.tree_data) {
      const int type = handler->type ();
      /* Pack type of tree data */
      mpiret = sc_MPI_Pack (&type, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
      SC_CHECK_MPI (mpiret);
      /* Pack each data */
      handler->pack_vector_prefix (buffer, num_bytes, pos, comm);
    }
  }

  void
  unpack (const void *buffer, const int num_bytes, int &pos, pseudo_tree &data, sc_MPI_Comm comm)
  {
    /* Clear existing tree data */
    data.tree_data.clear ();

    /* Unpack number of topological data */
    int topo_data_size = 0;
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &topo_data_size, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    data.topo_data.resize (topo_data_size);
    mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, data.topo_data.data (), topo_data_size, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    /* Unpack number of tree-specific data */
    int num_handler = 0;
    mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &num_handler, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    data.tree_data.resize (num_handler);
    for (auto &ihandler : data.tree_data) {
      int type;
      mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &type, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);

      if (type == 0) {
        ihandler = std::make_shared<t8_data_handler<enlarged_data<int>>> ();
      }
      else if (type == 1) {
        ihandler = std::make_shared<t8_data_handler<enlarged_data<double>>> ();
      }
      else {
        SC_ABORT_NOT_REACHED ();
      }

      int outcount = 0;
      ihandler->unpack_vector_prefix (buffer, num_bytes, pos, outcount, comm);
    }
  }

  int
  type ()
  {
    return -1;
  }
};

#endif /* T8_PSEUDO_TREES_HXX */
