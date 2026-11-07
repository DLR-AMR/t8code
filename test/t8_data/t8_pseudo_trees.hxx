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

/** \file t8_pseudo_trees.hxx
 * This file provides a pseudo tree structure for testing purposes.
 */

#ifndef T8_PSEUDO_TREES_HXX
#define T8_PSEUDO_TREES_HXX

#include <t8.h>
#include <vector>
#include <t8_data/t8_data_handler_specs.hxx>
#include <t8_data/t8_enlarged_stdtypes.hxx>
#include <t8_data/t8_vector_handler.hxx>
#include <memory>

/**
 * Represents a pseudo tree structure containing topological and tree data.
 *
 * \note This class is used for testing purposes only.
 */
struct pseudo_tree
{
  std::vector<int> topo_data; /**< A vector of integers representing topological data. */
  std::vector<std::shared_ptr<t8_abstract_vector_handler>>
    tree_data; /**< A vector of pointers to t8_abstract_vector_handler objects, representing tree data. */
};

/** Data handler for for the \ref pseudo_tree structure. */
template <>
struct t8_data_handler<pseudo_tree>
{
  /**
   * Returns the size of a pseudo_tree.
   *
   * \param[in] item The data to compute the size of.
   * \param[in] comm The MPI communicator used for communication.
   * \return An integer representing the size of the data.
   */
  int
  size (const pseudo_tree &item, sc_MPI_Comm comm) const
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

  /**
   * Packs a pseudo_tree into a buffer for communication.
   *
   * \param[in] item The data to be packed.
   * \param[in, out] buffer The buffer where the data will be packed.
   * \param[in] buffer_size The number of bytes available in the buffer.
   * \param[in, out] position The current position in the buffer. Advanced by the packed bytes.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  pack (const pseudo_tree &item, void *buffer, const int buffer_size, int &position, sc_MPI_Comm comm) const
  {
    const int data_size = item.topo_data.size ();
    /* Pack number of topological data */
    int mpiret = sc_MPI_Pack (&data_size, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack topological data in one call */
    mpiret = sc_MPI_Pack (item.topo_data.data (), data_size, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);
    /* Pack number of tree-specific data */
    const int tree_data_size = item.tree_data.size ();
    mpiret = sc_MPI_Pack (&tree_data_size, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
    SC_CHECK_MPI (mpiret);

    for (const auto &handler : item.tree_data) {
      const t8_data_handler_type type = handler->type ();
      /* Pack type of tree data */
      mpiret = sc_MPI_Pack (&type, 1, sc_MPI_INT, buffer, buffer_size, &position, comm);
      SC_CHECK_MPI (mpiret);
      /* Pack each data */
      handler->pack_vector_prefix (buffer, buffer_size, position, comm);
    }
  }

  /**
   * Unpacks a pseudo_tree from a buffer.
   *
   * \param[in] buffer A pointer to the buffer containing the packed data.
   * \param[in] buffer_size The number of bytes in the buffer.
   * \param[in, out] position The current position in the buffer. Advanced by the unpacked bytes.
   * \param[out] item The unpacked data.
   * \param[in] comm The MPI communicator used for communication.
   */
  void
  unpack (const void *buffer, const int buffer_size, int &position, pseudo_tree &item, sc_MPI_Comm comm) const
  {
    /* Clear existing tree data */
    item.tree_data.clear ();

    /* Unpack number of topological data */
    int topo_data_size = 0;
    int mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &topo_data_size, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);
    item.topo_data.resize (topo_data_size);
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, item.topo_data.data (), topo_data_size, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    /* Unpack number of tree-specific data */
    int num_handler = 0;
    mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &num_handler, 1, sc_MPI_INT, comm);
    SC_CHECK_MPI (mpiret);

    item.tree_data.resize (num_handler);
    for (auto &ihandler : item.tree_data) {
      t8_data_handler_type type;
      mpiret = sc_MPI_Unpack (buffer, buffer_size, &position, &type, 1, sc_MPI_INT, comm);
      SC_CHECK_MPI (mpiret);

      if (!t8_is_internal_data (type)) {

        /* TODO: This is currently only a placeholder for actual internal data types. */
        if (type.get () == 0) {
          ihandler = std::make_shared<t8_vector_handler<enlarged_data<int>>> ();
        }
        else if (type.get () == 1) {
          ihandler = std::make_shared<t8_vector_handler<enlarged_data<double>>> ();
        }
        else {
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {
        ihandler = t8_create_internal_handler (type);
      }
      int outcount = 0;
      ihandler->unpack_vector_prefix (buffer, buffer_size, position, outcount, comm);
    }
  }

  /**
   * Returns the type of the data handler.
   *
   * \return The type of the data handler.
   */
  t8_data_handler_type
  type () const
  {
    return t8_data_handler_type (-1);
  }
};

#endif /* T8_PSEUDO_TREES_HXX */
