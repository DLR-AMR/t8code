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

/** \file t8_mpi.h
 *
 *This file provides a consistent MPI interface for MPI datatypes
 * whether t8code is configure with MPI or without. This files
 * contains patched functions from libsc.
 */

#ifndef T8_MPI_H
#define T8_MPI_H

/* include config headers */
#ifndef T8_CMAKE_BUILD
#include <t8_config.h>
#endif

T8_EXTERN_C_BEGIN ();

#ifdef T8_ENABLE_MPI
#include <mpi.h>

#define t8_MPI_INT8_T MPI_INT8_T
#define t8_MPI_UNSIGNED_LONG_LONG MPI_UNSIGNED_LONG_LONG

#define t8_MPI_Type_size MPI_Type_size

#else /* !T8_ENABLE_MPI */

#define t8_MPI_INT8_T (1)
#define t8_MPI_UNSIGNED_LONG_LONG (0x4c000409)

/** Return size of an MPI datatype.
 * \param [in] datatype     Valid MPI datatype.
 * \param [out] size        Size of MPI datatype in bytes.
 * \return                  MPI_SUCCESS on success.
 */
int
t8_MPI_Type_size (sc_MPI_Datatype datatype, int *size);

#endif

/** Return the size of MPI datatypes.
 * \param [in] t    MPI datatype.
 * \return          Returns the size in bytes.
 */
size_t
t8_mpi_sizeof (sc_MPI_Datatype t);

T8_EXTERN_C_END ();

#endif /* !T8_MPI_H */
