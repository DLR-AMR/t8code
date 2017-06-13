/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_cmesh_zoltan.h
 *
 * We define routines to interface with the Zoltan library for mesh partitioning.
 * This header file should only be included, if t8code was configure with the
 * --with-zoltan option.
 *
 * TODO: document this file
 */

#ifndef T8_CMESH_ZOLTAN_H
#define T8_CMESH_ZOLTAN_H

#include <t8.h>

/* The function in this file can only be used if t8code was configured
 * with --with-zoltan.
 */
#ifdef T8_WITH_ZOLTAN
#include <zoltan.h>

#define T8_CHECK_ZOLTAN(r) SC_CHECK_ABORT ((r) == ZOLTAN_OK, "MPI error")

T8_EXTERN_C_BEGIN ();

/** Initialize t8code for use with Zoltan.
 * \param [in]  argc The number of command line arguments.
 * \param [in]  argc Array of command line arguments.
 * \note This function must be called once before working with Zoltan
 * routines.
 */
void                t8_cmesh_zoltan_initialize (int argc, char **argv);

/** Setup Zoltan for use with parmetis Graph partitioning methods for
 * a particular cmesh and store this settings at the cmesh.
 * \param [in,out] cmesh A committed cmesh.
 * \param [in]  comm  The MPI communicator that should be used. Must
 *                    fulfill \ref t8_cmesh_comm_is_valid.
 */
void                t8_cmesh_zoltan_setup_parmetis (t8_cmesh_t cmesh,
                                                    sc_MPI_Comm comm);

/** Free the memory used by the Zoltan methods.
 * \param [in] cmesh  A cmesh.
 * \note  \a cmesh should have been setup for Zoltan usage with t8_cmesh_zolten_setup_parmetis.
 */
void                t8_cmesh_zoltan_destroy (t8_cmesh_t cmesh);

T8_EXTERN_C_END ();

#endif

#endif /* !T8_CMESH_ZOLTAN_H */
