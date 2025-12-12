/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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
 * \file This file collects all functionality required to send elements from
 * one process to another within the partition-for-coarsening offset correction.
*/

#ifndef T8_FOREST_PFC_MESSAGE_H
#define T8_FOREST_PFC_MESSAGE_H

#include <t8.h>
#include <t8_element.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_pfc_helper.hxx>
#include <t8_schemes/t8_scheme.h>
#include <t8_schemes/t8_scheme.hxx>
#include <sc_mpi.h>

#define T8_PFC_MESSAGE 54323

/**
 * This class collects all functionality to handle the data to be sent/received
 * between the processes that will be required to decide on whether and where
 * families are split at process boundaries.
*/
class t8_forest_pfc_message {
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
  pack (void *buf, int buf_size, int *position)
  {
    /* pack: itree */
    sc_MPI_Pack (&itree, 1, T8_MPI_GLOIDX, buf, buf_size, position, comm);

    /* pack: eclass */
    int eclass_int = (int) eclass;
    sc_MPI_Pack (&eclass_int, 1, sc_MPI_INT, buf, buf_size, position, comm);

    /* pack: parent */
    t8_element_MPI_Pack (myscheme, eclass, &parent, 1, buf, buf_size, position, comm);

    /* pack: num_siblings */
    sc_MPI_Pack (&num_siblings, 1, sc_MPI_INT, buf, buf_size, position, comm);
  }

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
  unpack (void *buf, int buf_size, int *position)
  {
    /* unpack: itree */
    sc_MPI_Unpack (buf, buf_size, position, &itree, 1, T8_MPI_GLOIDX, comm);

    /* unpack: eclass */
    int eclass_int;
    sc_MPI_Unpack (buf, buf_size, position, &eclass_int, 1, sc_MPI_INT, comm);
    eclass = (t8_eclass_t) eclass_int;

    /* unpack: parent */
    t8_element_new (myscheme, eclass, 1, &parent);
    allocated_parent = true;
    t8_element_MPI_Unpack (myscheme, eclass, buf, buf_size, position, &parent, 1, comm);

    /* unpack: num_siblings */
    sc_MPI_Unpack (buf, buf_size, position, &num_siblings, 1, sc_MPI_INT, comm);
  }

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
  pack_size ()
  {
    /* initialize sum*/
    int message_size = 0;
    int datasize;

    /* add size: itree */
    sc_MPI_Pack_size (1, T8_MPI_GLOIDX, comm, &datasize);
    message_size += datasize;

    /* add size: eclass */
    sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    message_size += datasize;

    /* add size: parent */
    if (parent != NULL) {
      // t8_eclass_scheme_c *scheme = schemes->eclass_schemes[eclass];
      // scheme->t8_element_MPI_Pack_size (1, comm, &datasize);
      t8_element_MPI_Pack_size (myscheme, eclass, 1, comm, &datasize);
      message_size += datasize;
    }

    /* add size: num_siblings */
    sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    message_size += datasize;
    return message_size;
  }

  /**
   * Send the message using non-blocking MPI send.
   *
   * \param[in]   forest  the forest (only used for forest->mpicomm)
   * \param[out]  request the MPI send request
  */
  void
  mpi_Isend (const t8_forest_t forest, sc_MPI_Request &request)
  {
    /* Allocate buffer */
    int buffer_size = pack_size ();
    char *send_buffer = T8_ALLOC (char, buffer_size);

    /* Pack message to buffer */
    int position = 0;
    pack (send_buffer, buffer_size, &position);

    /* Send buffer */
    const int mpiret
      = sc_MPI_Isend (send_buffer, position, sc_MPI_PACKED, iproc, message_tag, forest->mpicomm, &request);
    SC_CHECK_MPI (mpiret);
    T8_FREE (send_buffer);
  }

  /**
   * Receive data from another process. Probes for the size of the message before receiving
   *
   * \param[in, out]  recv_buf  the receive buffer. The function allocates the right size of memory to receive the message. 
   * \param[out]  buf_size  the size of the buffer
  */
  void
  mpi_Recv (char *&recv_buf, int &buf_size)
  {
    /** Get needed size of message via MPI probe and allocate buffer. */
    sc_MPI_Status status;
    sc_MPI_Probe (iproc, message_tag, comm, &status);
    sc_MPI_Get_count (&status, sc_MPI_PACKED, &buf_size);
    recv_buf = T8_ALLOC (char, buf_size);

    /* Actually receive buffer. */
    const int mpiret = sc_MPI_Recv (recv_buf, buf_size, sc_MPI_PACKED, iproc, message_tag, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /**
   * Fill this instance of class t8_forest_pfc_message, i.e., setting the member viarables based on the given forest.
   *
   * \param[in] forest the forest
  */
  void
  fill (t8_forest_t forest)
  {
    // Set mpi rank and partition element offsets.
    const t8_procidx_t rank = forest->mpirank;
    const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);

    // Determine the global element ID of the current process ("rank") that is closest to the elements of process iproc:
    // If iproc < rank, the first element is closest to iproc; otherwise, the last one is.
    const t8_gloidx_t closest_to_rank_gid = (iproc <= rank) ? partition[rank] : partition[rank + 1] - 1;

    // Use helper function to get various element and tree indices.
    t8_tree_t tree;
    t8_locidx_t index_in_tree;
    t8_element_t *element_closest_to_receiver;
    t8_forest_pfc_helper_index_in_tree_from_globalid (forest, closest_to_rank_gid, itree, tree, index_in_tree,
                                                      element_closest_to_receiver);
    // Set scheme and eclass.
    myscheme = t8_forest_get_scheme (forest);
    eclass = t8_forest_get_eclass (forest, t8_forest_get_local_id (forest, itree));

    // If we are already the root element, we cannot be part of a split family, so we send any(the root) element and no num_siblings.
    if (myscheme->element_get_level (eclass, element_closest_to_receiver) == 0) {
      parent = element_closest_to_receiver;
      num_siblings = 0;
    }
    else {
      // Compute parent.
      t8_element_new (myscheme, eclass, 1, &parent);
      allocated_parent = true;
      myscheme->element_get_parent (eclass, element_closest_to_receiver, parent);

      // Distinguish send "direction"
      if (iproc > rank) {
        // Sending will be towards higher-rank process, so count siblings in decreasing index direction.
        const t8_locidx_t min_id = t8_forest_pfc_extreme_local_sibling (myscheme, tree, index_in_tree, true);
        num_siblings = index_in_tree - min_id + 1;
      }
      else {
        T8_ASSERT (iproc != rank);
        // Sending will be towards lower-rank process, so count siblings in increasing index direction.
        const t8_locidx_t max_id = t8_forest_pfc_extreme_local_sibling (myscheme, tree, index_in_tree, false);
        num_siblings = max_id - index_in_tree + 1;
      }
    }
  }

  /**
   * Getter function for the parent to stay in control of memory management.
   *
   * \return The parent to be sent.
  */
  const t8_element_t *
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
    : itree (0), eclass (T8_ECLASS_ZERO), num_siblings (0), myscheme (scheme), comm (comm), iproc (iproc),
      message_tag (T8_PFC_MESSAGE), parent (NULL), allocated_parent (0)
  {
  }

  /// No (implicit) copy constructor.
  t8_forest_pfc_message (const t8_forest_pfc_message &other) = delete;

  /**
   * Move constructor.
  */
  t8_forest_pfc_message (t8_forest_pfc_message &&other)
    : itree { other.itree }, eclass (other.eclass), num_siblings (other.num_siblings), myscheme (other.myscheme),
      comm (other.comm), iproc (other.iproc), parent (other.parent), allocated_parent (other.allocated_parent)
  {
    if (allocated_parent) {
      other.parent = NULL;
      other.allocated_parent = false;
    }
  }

  /**
   * Default destructor
  */
  ~t8_forest_pfc_message ()
  {
    if (allocated_parent) {
      myscheme->element_destroy (eclass, 1, &parent);
      parent = NULL; /*TODO: necessary?*/
    }
  }

  // Data directly sent in the message:
  t8_gloidx_t itree;  /**< the global tree id */
  t8_eclass_t eclass; /**< the tree's eclass */
  int num_siblings;   /**< the process-local number of siblings */

  // Auxiliary data:
  const t8_scheme_c *myscheme; /**< the scheme class */
  sc_MPI_Comm comm;            /**< the MPI communicator */
  t8_procidx_t iproc;          /**< the process to send data to */
  int message_tag;             /**< the TAG identifying the message */

 private:
  t8_element_t *parent;  /**< The parent element to be sent. */
  bool allocated_parent; /**< Are we responsible for memory management of parent? */
};

#endif
