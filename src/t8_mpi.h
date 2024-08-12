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

#define t8_MPI_Pack MPI_Pack
#define t8_MPI_Pack_size MPI_Pack_size
#define t8_MPI_Unpack MPI_Unpack

#define t8_MPI_Gather MPI_Gather
#define t8_MPI_Allgather MPI_Allgather

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

/** Determine space needed to pack several instances of the same datatype.
 * \param [in] incount        Number of elements to pack.
 * \param [in] datatype       Datatype of elements to pack.
 * \param [in] comm           Valid MPI communicator.
 * \param [out] size          Number of bytes needed to packed \b incount
 *                            instances of \b datatype.
 * \return                    MPI_SUCCESS on success.
 */
int
t8_MPI_Pack_size (int incount, sc_MPI_Datatype datatype, sc_MPI_Comm comm, int *size);

/** Pack several instances of the same datatype into contiguous memory.
 * \param [in] inbuf          Buffer of elements of type \b datatype.
 * \param [in] incount        Number of elements in \b inbuf.
 * \param [in] datatype       Datatype of elements in \b inbuf.
 * \param [out] outbuf        Output buffer in which elements are packed.
 * \param [in] outsize        Size of output buffer in bytes.
 * \param [in, out] position  The current position in the output buffer.
 * \param [in] comm           Valid MPI communicator.
 * \return                    MPI_SUCCESS on success.
 */
int
t8_MPI_Pack (const void *inbuf, int incount, sc_MPI_Datatype datatype, void *outbuf, int outsize, int *position,
             sc_MPI_Comm comm);

/** Unpack contiguous memory into several instances of the same datatype.
 * \param [in] inbuf          Buffer of packed data.
 * \param [in] insize         Number of bytes in \b inbuf
 * \param [in, out] position  The current position in the input buffer.
 * \param [out] outbuf        Output buffer in which elements are unpacked.
 * \param [in] outcount       Number of elements to unpack.
 * \param [in] datatype       Datatype of elements to be unpacked.
 * \param [in] comm           Valid MPI communicator.
 * \return                    MPI_SUCCESS on success.
 */
int
t8_MPI_Unpack (const void *inbuf, int insize, int *position, void *outbuf, int outcount, sc_MPI_Datatype datatype,
               sc_MPI_Comm comm);

int
t8_MPI_Gather (void *, int, sc_MPI_Datatype, void *, int, sc_MPI_Datatype, int, sc_MPI_Comm);

/** Execute the MPI_Allgather algorithm. */
int
sc_MPI_Allgather (void *, int, sc_MPI_Datatype, void *, int, sc_MPI_Datatype, sc_MPI_Comm);
#endif

/** Return the size of MPI datatypes.
 * \param [in] t    MPI datatype.
 * \return          Returns the size in bytes.
 */
size_t
t8_mpi_sizeof (sc_MPI_Datatype t);

/** Fill a shmem array with an allgather.
 *
 * \param[in] sendbuf         the source from this process
 * \param[in] sendcount       the number of items to allgather
 * \param[in] sendtype        the type of items to allgather
 * \param[in,out] recvbuf     the destination shmem array
 * \param[in] recvcount       the number of items to allgather
 * \param[in] recvtype        the type of items to allgather
 * \param[in] comm            the mpi communicator
 */
void
t8_shmem_allgather (void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, void *recvbuf, int recvcount,
                    sc_MPI_Datatype recvtype, sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* !T8_MPI_H */
