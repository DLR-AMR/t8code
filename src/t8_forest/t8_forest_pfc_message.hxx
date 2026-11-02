/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/**
 * \file t8_forest_pfc_message.hxx
 * This file collects all functionality required to send elements from
 * one process to another within the partition-for-coarsening offset correction.
*/

#ifndef T8_FOREST_PFC_MESSAGE_H
#define T8_FOREST_PFC_MESSAGE_H

#include <t8.h>
#include <t8_element/t8_element.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <sc_mpi.h>

/**
 * This class collects all functionality to handle the data to be sent/received
 * between the processes that will be required to decide on whether and where
 * families are split at process boundaries.
*/
struct t8_forest_pfc_message
{
 public:
  /**
   * Pack the data to prepare sending.
   *
   * \param[in,out] buf       the sending buffer
   *                            on input: allocated (but empty)
   *                            on output: filled with the data to sent
   * \param[in]     buf_size  the size of the sending buffer
   * \param[in,out] position  current position within the sending buffer
   *                            on input:  0
   *                            on output: buf_size
   */
  void
  pack (void *buf, int buf_size, int *position);

  /**
   * Unpack the received data.
   *
   * Transfer the received data into the member variables of this class.
   *
   * \param[in]     buf       the buffer containing the received MPI data
   * \param[in]     buf_size  the size of the buffer
   * \param[in,out] position  current position within the receive buffer
   *                            on input:  0
   *                            on output: buf_size
   */
  void
  unpack (void *buf, int buf_size, int *position);

  /**
   * Determine and return the pack size of the PFC message.
   *
   * It results from the data fields to send and the associated sizes:
   *  + itree (T8_MPI_GLOIDX)
   *  + eclass (sc_MPI_INT)
   *  + parent (via t8_element_MPI_Pack_size)
   *  + num_siblings (sc_MPI_INT)
   *
   * \return The pack size of the message.
  */
  int
  pack_size ();

  /**
   * Send the message using non-blocking MPI send.
   *
   * \param[in]   forest  the forest (only used for forest->mpicomm)
   * \param[out]  request the MPI send request
  */
  void
  mpi_Isend (const t8_forest_t forest, sc_MPI_Request &request);

  /**
   * Receive data from another process. Probes for the size of the message before receiving
   *
   * \param[in, out]  recv_buf  the receive buffer. The function allocates the right size of memory to receive the message.
   * \param[out]  buf_size  the size of the buffer
  */
  void
  mpi_Recv (char *&recv_buf, int &buf_size);

  /**
   * Fill this instance of class t8_forest_pfc_message, i.e., setting the member viarables based on the given forest.
   *
   * \param[in] forest the forest
  */
  void
  fill (t8_forest_t forest);

  /**
   * Getter function for the parent to stay in control of memory management.
   *
   * \return The parent to be sent.
  */
  inline const t8_element_t *
  get_parent () const
  {
    T8_ASSERT (parent);
    return parent;
  }

  /**
   * Constructor of class t8_forest_pfc_message.
   *
   * The arguments are directly copied into the corresponding member variables;
   * the remaining member variables are set to default values.
   *
   * \param[in] scheme the scheme
   * \param[in] iproc     the process to send to
   * \param[in] comm      the MPI communicator
   *
  */
  t8_forest_pfc_message (const t8_scheme_c *scheme, t8_procidx_t iproc, sc_MPI_Comm comm)
    : itree (0), eclass (T8_ECLASS_ZERO), num_siblings (0), scheme (scheme), comm (comm), iproc (iproc),
      parent (nullptr)
  {
  }

  /// No (implicit) copy constructor.
  t8_forest_pfc_message (const t8_forest_pfc_message &other) = delete;

  /**
   * Move constructor.
  */
  t8_forest_pfc_message (t8_forest_pfc_message &&other)
    : itree { other.itree }, eclass (other.eclass), num_siblings (other.num_siblings), scheme (other.scheme),
      comm (other.comm), iproc (other.iproc), parent (other.parent)
  {
    if (parent != nullptr) {
      other.parent = nullptr;
    }
  }

  /**
   * Default destructor
  */
  ~t8_forest_pfc_message ()
  {
    if (parent != nullptr) {
      scheme->element_destroy (eclass, 1, &parent);
      parent = nullptr; /*TODO: necessary?*/
    }
  }

  // Data directly sent in the message:
  t8_gloidx_t itree;  /**< the global tree id */
  t8_eclass_t eclass; /**< the tree's eclass */
  int num_siblings;   /**< the process-local number of siblings */

  // Auxiliary data:
  const t8_scheme_c *scheme; /**< the scheme class */
  sc_MPI_Comm comm;          /**< the MPI communicator */
  t8_procidx_t iproc;        /**< the process to send data to */

 private:
  t8_element_t *parent; /**< The parent element to be sent. */
};

#endif
