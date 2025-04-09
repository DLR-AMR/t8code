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
 * \file t8_cmesh_examples.h
 * We provide example coarse meshes in this file 
 */

#ifndef T8_CMESH_EXAMPLES
#define T8_CMESH_EXAMPLES
#include <t8_cmesh.h>
#include <p4est_connectivity.h>
#include <p8est_connectivity.h>
#include <t8_cmesh/t8_cmesh_geometry.h>

T8_EXTERN_C_BEGIN ();

/** Constructs a cmesh from a given p4est_connectivity structure.
 * \param[in]       conn       The p4est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 * \note This function requires that p4est is initialized. Make sure to call
 * \ref p4est_init before using this routine. If this is not the case, a
 * warning is issued and \ref p4est_init is called from within this function.
 */
t8_cmesh_t
t8_cmesh_new_from_p4est (p4est_connectivity_t *conn, sc_MPI_Comm comm, int do_partition);

/** Constructs a cmesh from a given p8est_connectivity structure.
 * \param[in]       conn       The p8est connectivity.
 * \param[in]       comm       mpi communicator to be used with the new cmesh.
 * \param[in]       do_dup     Flag whether the communicator shall be duplicated or not.
 * \param[in]       do_partition Flag whether the cmesh should be partitioned or not.
 * \return          A t8_cmesh structure that holds the same connectivity information
 *                  as \a conn.
 * \note This function requires that p4est is initialized. Make sure to call
 * \ref p4est_init before using this routine. If this is not the case, a
 * warning is issued and \ref p4est_init is called from within this function.
 */
t8_cmesh_t
t8_cmesh_new_from_p8est (p8est_connectivity_t *conn, sc_MPI_Comm comm, int do_partition);

T8_EXTERN_C_END ();

#endif /* !T8_CMESH_EXAMPLES */
