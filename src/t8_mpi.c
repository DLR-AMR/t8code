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

#include <t8.h>
#include <t8_mpi.h>

size_t
t8_mpi_sizeof (sc_MPI_Datatype t)
{
  if (t == sc_MPI_BYTE || t == t8_MPI_INT8_T)
    return 1;
  if (t == sc_MPI_CHAR || t == sc_MPI_UNSIGNED_CHAR)
    return sizeof (char);
  if (t == sc_MPI_SHORT || t == sc_MPI_UNSIGNED_SHORT)
    return sizeof (short);
  if (t == sc_MPI_INT || t == sc_MPI_UNSIGNED)
    return sizeof (int);
  if (t == sc_MPI_LONG || t == sc_MPI_UNSIGNED_LONG)
    return sizeof (long);
  if (t == sc_MPI_LONG_LONG_INT || t == t8_MPI_UNSIGNED_LONG_LONG)
    return sizeof (long long);
  if (t == sc_MPI_FLOAT)
    return sizeof (float);
  if (t == sc_MPI_DOUBLE)
    return sizeof (double);
  if (t == sc_MPI_LONG_DOUBLE)
    return sizeof (long double);
  if (t == sc_MPI_2INT)
    return 2 * sizeof (int);

  SC_ABORT_NOT_REACHED ();
}

#ifndef T8_ENABLE_MPI

int
t8_MPI_Type_size (sc_MPI_Datatype datatype, int *size)
{
  *size = t8_mpi_sizeof (datatype);

  return sc_MPI_SUCCESS;
}

int
t8_MPI_Pack_size (int incount, sc_MPI_Datatype datatype, sc_MPI_Comm comm, int *size)
{
  int mpiret;

  SC_ASSERT (incount >= 0);
  SC_ASSERT (size != NULL);

  mpiret = t8_MPI_Type_size (datatype, size);
  SC_CHECK_MPI (mpiret);
  *size *= incount;

  return sc_MPI_SUCCESS;
}

int
t8_MPI_Pack (const void *inbuf, int incount, sc_MPI_Datatype datatype, void *outbuf, int outsize, int *position,
             sc_MPI_Comm comm)
{
  int mpiret;
  int size;

  SC_ASSERT (incount >= 0);
  SC_ASSERT (position != NULL);

  mpiret = t8_MPI_Pack_size (incount, datatype, comm, &size);
  SC_CHECK_MPI (mpiret);

  /* Check that we have enough space to pack the datatypes */
  if (*position + size > outsize) {
    return sc_MPI_ERR_NO_SPACE;
  }

  /* Copy the contiguous memory */
  memcpy ((char *) outbuf + *position, inbuf, size);
  *position += size;

  return sc_MPI_SUCCESS;
}

int
t8_MPI_Unpack (const void *inbuf, int insize, int *position, void *outbuf, int outcount, sc_MPI_Datatype datatype,
               sc_MPI_Comm comm)
{
  int mpiret;
  int size;

  SC_ASSERT (position != NULL);
  SC_ASSERT (outcount >= 0);

  mpiret = t8_MPI_Pack_size (outcount, datatype, comm, &size);
  SC_CHECK_MPI (mpiret);

  /* Check that the message is big enough for the datatypes that we want */
  if (*position + size > insize) {
    return sc_MPI_ERR_NO_SPACE;
  }

  /* Copy the contiguous memory */
  memcpy (outbuf, (char *) inbuf + *position, size);
  *position += size;

  return sc_MPI_SUCCESS;
}

int
t8_MPI_Gather (void *p, int np, sc_MPI_Datatype tp, void *q, int nq, sc_MPI_Datatype tq, int rank, sc_MPI_Comm comm)
{
  size_t lp;
#ifdef SC_ENABLE_DEBUG
  size_t lq;
#endif

  SC_ASSERT (rank == 0 && np >= 0 && nq >= 0);

  /* *INDENT-OFF* horrible indent bug */
  lp = (size_t) np * t8_mpi_sizeof (tp);
#ifdef SC_ENABLE_DEBUG
  lq = (size_t) nq * t8_mpi_sizeof (tq);
#endif
  /* *INDENT-ON* */

  SC_ASSERT (lp == lq);
  memcpy (q, p, lp);

  return sc_MPI_SUCCESS;
}

#endif
