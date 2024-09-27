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

#include <memory>

/**
 * \class pseudo_tree
 * A class representing a pseudo tree structure.
 *
 * The pseudo_tree class encapsulates a tree-like structure with topology data and tree data.
 * It provides constructors, assignment operator, and destructor for managing the tree data.
 * 
 * It is meant for testing purposes only and to mimic the structure of a tree.
 *
 * \var std::vector<int> pseudo_tree::topo_data
 * A vector containing the topology data of the pseudo tree.
 *
 * \var std::vector<std::unique_ptr<t8_abstract_data_handler>> pseudo_tree::tree_data
 * A vector of unique pointers to t8_abstract_data_handler objects representing the tree data.
 */
class pseudo_tree {
 public:
  pseudo_tree ()
  {
  }

  pseudo_tree (const pseudo_tree &other): topo_data (other.topo_data)
  {
    tree_data.resize (other.tree_data.size ());
    for (size_t i = 0; i < other.tree_data.size (); ++i) {
      tree_data[i] = std::unique_ptr<t8_abstract_data_handler> (other.tree_data[i]->clone ());
    }
  }

  pseudo_tree &
  operator= (const pseudo_tree &other)
  {
    if (this != &other) {
      tree_data.clear ();
      topo_data = other.topo_data;

      tree_data.resize (other.tree_data.size ());
      for (size_t i = 0; i < other.tree_data.size (); ++i) {
        tree_data[i] = std::unique_ptr<t8_abstract_data_handler> (other.tree_data[i]->clone ());
      }
    }
    return *this;
  }

  ~pseudo_tree () = default;

  std::vector<int> topo_data;
  std::vector<std::unique_ptr<t8_abstract_data_handler>> tree_data;
};

template <>
class t8_single_data_handler<pseudo_tree> {
 public:
  int
  size (const pseudo_tree *item, sc_MPI_Comm comm)
  {
    int int_size = 0;
    const int topo_data_size = item->topo_data.size () + 1;

    int mpiret = sc_MPI_Pack_size (1, sc_MPI_INT, comm, &int_size);
    SC_CHECK_MPI (mpiret);
    int total_size = topo_data_size * int_size;

    /* tree_data_size */
    total_size += int_size;
    for (const auto &ihandler : item->tree_data) {
      total_size += ihandler->buffer_size (comm) + int_size;
    }
    return total_size;
  }

  void
  pack (const pseudo_tree *data, int &pos, void *buffer, const int num_bytes, sc_MPI_Comm comm)
  {
    const int data_size = data->topo_data.size ();
    /* Pack number of topological data */
    int mpiret = sc_MPI_Pack (&data_size, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack each topological data*/
    mpiret = sc_MPI_Pack ((data->topo_data.data ()), data_size, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack number of tree-specific data*/
    const int tree_data_size = data->tree_data.size ();
    mpiret = sc_MPI_Pack (&tree_data_size, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
    SC_CHECK_MPI (mpiret);

    for (auto &handler : data->tree_data) {
      const int type = handler->type ();
      /* Pack type of tree data */
      mpiret = sc_MPI_Pack (&type, 1, sc_MPI_INT, buffer, num_bytes, &pos, comm);
      SC_CHECK_MPI (mpiret);
      /* Pack each data. */
      handler->pack_vector_prefix (buffer, num_bytes, pos, comm);
    }
  }

  void
  unpack (const void *buffer, const int num_bytes, int &pos, pseudo_tree *data, sc_MPI_Comm comm)
  {
    /* Clear existing tree data */
    for (const auto &handler_ptr : data->tree_data) {
      t8_abstract_data_handler *handler = handler_ptr.get ();
      delete handler;
    }
    data->tree_data.clear ();

    /* Unpack number of topological data */
    int topo_data_size = 0;
    int mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &topo_data_size, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    data->topo_data.resize (topo_data_size);
    for (int &topo_item : data->topo_data) {
      /* Unpack each topological item */
      mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &topo_item, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);
    }
    /* Unpack number of tree-specific data */
    int num_handler = 0;
    mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &num_handler, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    for (int ihandler = 0; ihandler < num_handler; ihandler++) {
      int type;
      mpiret = sc_MPI_Unpack (buffer, num_bytes, &pos, &type, 1, sc_MPI_INT, comm);
      int outcount = 0;
      if (type == 0) {
        t8_data_handler<enlarged_data<int>> *new_handler = new t8_data_handler<enlarged_data<int>> ();
        new_handler->unpack_vector_prefix (buffer, num_bytes, pos, outcount, comm);
        data->tree_data.push_back (std::unique_ptr<t8_abstract_data_handler> (new_handler));
      }
      else if (type == 1) {
        t8_data_handler<enlarged_data<double>> *new_handler = new t8_data_handler<enlarged_data<double>> ();
        new_handler->unpack_vector_prefix (buffer, num_bytes, pos, outcount, comm);
        data->tree_data.push_back (std::unique_ptr<t8_abstract_data_handler> (new_handler));
      }
      else {
        SC_ABORT_NOT_REACHED ();
      }
    }
  }

  int
  type ()
  {
    return -1;
  }
};

#endif /* T8_PSEUDO_TREES_HXX */
