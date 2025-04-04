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
class t8_forest_pfc_message_c {
 public:
  void
  pack (void *buf, int buf_size, int *position)
  {
    /* itree */
    sc_MPI_Pack (&itree, 1, T8_MPI_GLOIDX, buf, buf_size, position, comm);

    /* eclass */
    int eclass_int = (int) eclass;
    sc_MPI_Pack (&eclass_int, 1, sc_MPI_INT, buf, buf_size, position, comm);

    /* parent */
    // t8_eclass_scheme_c *scheme = schemes->eclass_schemes[eclass_int];
    // scheme->t8_element_MPI_Pack (parent, 1, buf, buf_size, position, comm);
    t8_element_MPI_Pack (myscheme, eclass, &parent, 1, buf, buf_size, position, comm);

    sc_MPI_Pack (&num_siblings, 1, sc_MPI_INT, buf, buf_size, position, comm);
  }

  void
  unpack (void *buf, int buf_size, int *position)
  {
    /* itree */
    sc_MPI_Unpack (buf, buf_size, position, &itree, 1, T8_MPI_GLOIDX, comm);

    /* eclass */
    int eclass_int;
    sc_MPI_Unpack (buf, buf_size, position, &eclass_int, 1, sc_MPI_INT, comm);
    eclass = (t8_eclass_t) eclass_int;

    /* parent */
    // t8_eclass_scheme_c *scheme = schemes->eclass_schemes[eclass_int];
    // scheme->t8_element_new (1, &parent);
    t8_element_new (myscheme, eclass, 1, &parent);
    allocated_parent = true;
    // scheme->t8_element_MPI_Unpack (buf, buf_size, position, parent, 1, comm);
    t8_element_MPI_Unpack (myscheme, eclass, buf, buf_size, position, &parent, 1, comm);

    /* num_sibling */
    sc_MPI_Unpack (buf, buf_size, position, &num_siblings, 1, sc_MPI_INT, comm);
  }

  int
  pack_size ()
  {
    int message_size = 0;
    int datasize;

    /* itree */
    sc_MPI_Pack_size (1, T8_MPI_GLOIDX, comm, &datasize);
    message_size += datasize;

    /* eclass */
    sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    message_size += datasize;

    /* parent */
    if (parent != NULL) {
      // t8_eclass_scheme_c *scheme = schemes->eclass_schemes[eclass];
      // scheme->t8_element_MPI_Pack_size (1, comm, &datasize);
      t8_element_MPI_Pack_size (myscheme, eclass, 1, comm, &datasize);
      message_size += datasize;
    }

    /* num_sibling */
    sc_MPI_Pack_size (1, sc_MPI_INT, comm, &datasize);
    message_size += datasize;
    return message_size;
  }

  void
  mpi_Isend (t8_forest_t forest, sc_MPI_Request &request)
  {
    /* Allocate buffer */
    int buffer_size = pack_size ();
    char *send_buffer = T8_ALLOC (char, buffer_size);

    /* Pack message to buffer */
    int position = 0;
    pack (send_buffer, buffer_size, &position);

    /* Send buffer */
    int mpiret = sc_MPI_Isend (send_buffer, position, sc_MPI_PACKED, iproc, message_tag, forest->mpicomm, &request);
    SC_CHECK_MPI (mpiret);
    T8_FREE (send_buffer);
  }

  void
  mpi_Recv (char *&recv_buf, int &buf_size)
  {
    /** Get needed size for packed message, and allocate */
    sc_MPI_Status status;
    sc_MPI_Probe (iproc, message_tag, comm, &status);
    sc_MPI_Get_count (&status, sc_MPI_PACKED, &buf_size);
    recv_buf = T8_ALLOC (char, buf_size);

    /* Actually receive buffer */
    int mpiret = sc_MPI_Recv (recv_buf, buf_size, sc_MPI_PACKED, iproc, message_tag, comm, sc_MPI_STATUS_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  void
  fill (t8_forest_t forest)
  {
    t8_procidx_t rank = forest->mpirank;
    const t8_gloidx_t *partition = t8_shmem_array_get_gloidx_array (forest->element_offsets);
    t8_gloidx_t closest_to_rank_gid = (iproc <= rank) ? partition[rank] : partition[rank + 1] - 1;
    t8_tree_t tree;
    t8_locidx_t index_in_tree;
    t8_element_t *element_closest_to_receiver;
    // t8_eclass_scheme_c *scheme;
    t8_forest_pfc_helper_index_in_tree_from_globalid (forest, closest_to_rank_gid, itree, tree, index_in_tree,
                                                      element_closest_to_receiver);

    // eclass = scheme->eclass;
    myscheme = t8_forest_get_scheme (forest);
    eclass = t8_forest_get_eclass (forest, itree);
    /* if we are already the root element, we cannot be part of a split family, so we send any(the root) element and no num_siblings */
    if (myscheme->element_get_level (eclass, element_closest_to_receiver) == 0) {
      parent = element_closest_to_receiver;
      num_siblings = 0;
    }
    else {
      /* compute parent and num_siblings*/
      // scheme->t8_element_new (1, &parent);
      t8_element_new (myscheme, eclass, 1, &parent);

      allocated_parent = true;
      // scheme->t8_element_parent (element_closest_to_receiver, parent);
      myscheme->element_get_parent (eclass, element_closest_to_receiver, parent);

      if (iproc > rank) {
        t8_locidx_t min_id = t8_forest_pfc_extreme_local_sibling (myscheme, tree, index_in_tree, true);
        num_siblings = index_in_tree - min_id + 1;
      }
      else {
        T8_ASSERT (iproc != rank);
        t8_locidx_t max_id = t8_forest_pfc_extreme_local_sibling (myscheme, tree, index_in_tree, false);
        num_siblings = max_id - index_in_tree + 1;
      }
    }
  }

  /* getter function to stay in control of memory management */
  const t8_element_t *
  get_parent ()
  {
    T8_ASSERT (parent);
    return parent;
  }
  // t8_forest_pfc_message_c (t8_scheme_cxx_t *schemes, t8_procidx_t iproc, sc_MPI_Comm comm)
  t8_forest_pfc_message_c (const t8_scheme_c *newscheme, t8_procidx_t iproc, sc_MPI_Comm comm)
    : itree (0), eclass (T8_ECLASS_ZERO), num_siblings (0), myscheme (newscheme), comm (comm), iproc (iproc),
      message_tag (T8_PFC_MESSAGE), parent (NULL), allocated_parent (0)
  {
    // t8_scheme_cxx_ref (myscheme);
    // t8_scheme_ref(myscheme);
  }
  t8_forest_pfc_message_c (const t8_forest_pfc_message_c &other) = delete;
  t8_forest_pfc_message_c (t8_forest_pfc_message_c &&other)
    : itree { other.itree }, eclass (other.eclass), num_siblings (other.num_siblings), myscheme (other.myscheme),
      comm (other.comm), iproc (other.iproc), parent (other.parent), allocated_parent (other.allocated_parent)
  {
    if (allocated_parent) {
      other.parent = NULL;
      other.allocated_parent = false;
    }
  }
  ~t8_forest_pfc_message_c ()
  {
    if (allocated_parent) {
      // t8_eclass_scheme_c *scheme = schemes->eclass_schemes[eclass];
      // scheme->t8_element_destroy (1, &parent);
      myscheme->element_destroy (eclass, 1, &parent);

      parent = NULL; /*necessary?*/
      // MYTODO
      // t8_scheme_cxx_unref (&schemes);
      // t8_scheme_unref(&myscheme);
    }
  }

  t8_gloidx_t itree;  /*directly sent*/
  t8_eclass_t eclass; /*directly sent*/
  int num_siblings;   /*directly sent*/

  /* background information -> to baseclass? */
  // t8_scheme_cxx_t *schemes;

  const t8_scheme_c *myscheme;
  // t8_eclass_t *my_tree_class;

  sc_MPI_Comm comm;
  t8_procidx_t iproc;
  int message_tag;

 private:
  t8_element_t *parent;  /*is unpacked and then sent*/
  bool allocated_parent; /* are we responsible for memory management of parent? */
};

#endif
