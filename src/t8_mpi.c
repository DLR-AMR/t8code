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

/* including t8_mpi.h does not work here since t8_mpi.h is included by t8.h */
#include <t8.h>
#include <sc_shmem.h>

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

int
t8_MPI_Allgather (void *p, int np, sc_MPI_Datatype tp, void *q, int nq, sc_MPI_Datatype tq, sc_MPI_Comm comm)
{
  return t8_MPI_Gather (p, np, tp, q, nq, tq, 0, comm);
}

#endif

static sc_shmem_type_t
t8_shmem_get_type_default (sc_MPI_Comm comm)
{
  sc_shmem_type_t type = sc_shmem_get_type (comm);
  if (type == SC_SHMEM_NOT_SET) {
    type = sc_shmem_default_type;
    sc_shmem_set_type (comm, type);
  }
  return type;
}

static void
t8_shmem_allgather_basic (void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, void *recvbuf, int recvcount,
                          sc_MPI_Datatype recvtype, sc_MPI_Comm comm, sc_MPI_Comm intranode, sc_MPI_Comm internode)
{
  int mpiret = t8_MPI_Allgather (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  SC_CHECK_MPI (mpiret);
}

#if defined(__bgq__) || defined(SC_ENABLE_MPIWINSHARED)
static void
t8_shmem_allgather_common (void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, void *recvbuf, int recvcount,
                           sc_MPI_Datatype recvtype, sc_MPI_Comm comm, sc_MPI_Comm intranode, sc_MPI_Comm internode)
{
  size_t typesize;
  int mpiret, intrarank, intrasize;
  char *noderecvchar = NULL;

  typesize = t8_mpi_sizeof (recvtype);

  mpiret = sc_MPI_Comm_rank (intranode, &intrarank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (intranode, &intrasize);
  SC_CHECK_MPI (mpiret);

  /* node root gathers from node */
  if (!intrarank) {
    noderecvchar = SC_ALLOC (char, intrasize *recvcount *typesize);
  }
  mpiret = t8_MPI_Gather (sendbuf, sendcount, sendtype, noderecvchar, recvcount, recvtype, 0, intranode);
  SC_CHECK_MPI (mpiret);

  /* node root allgathers between nodes */
  if (sc_shmem_write_start (recvbuf, comm)) {
    mpiret = t8_MPI_Allgather (noderecvchar, sendcount * intrasize, sendtype, recvbuf, recvcount * intrasize, recvtype,
                               internode);
    SC_CHECK_MPI (mpiret);
    SC_FREE (noderecvchar);
  }
  sc_shmem_write_end (recvbuf, comm);
}
#endif

#if !defined(SC_SHMEM_DEFAULT)
#define SC_SHMEM_DEFAULT SC_SHMEM_BASIC
#endif

void
t8_shmem_allgather (void *sendbuf, int sendcount, sc_MPI_Datatype sendtype, void *recvbuf, int recvcount,
                    sc_MPI_Datatype recvtype, sc_MPI_Comm comm)
{
  sc_shmem_type_t type;
  sc_MPI_Comm intranode = sc_MPI_COMM_NULL, internode = sc_MPI_COMM_NULL;

  type = t8_shmem_get_type_default (comm);
  sc_mpi_comm_get_node_comms (comm, &intranode, &internode);
  if (intranode == sc_MPI_COMM_NULL || internode == sc_MPI_COMM_NULL) {
    type = SC_SHMEM_BASIC;
  }
  switch (type) {
  case SC_SHMEM_BASIC:
  case SC_SHMEM_PRESCAN:
    t8_shmem_allgather_basic (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, intranode, internode);
    break;
#if defined(__bgq__) || defined(SC_ENABLE_MPIWINSHARED)
#if defined(__bgq__)
  case SC_SHMEM_BGQ:
  case SC_SHMEM_BGQ_PRESCAN:
#endif
#if defined(SC_ENABLE_MPIWINSHARED)
  case SC_SHMEM_WINDOW:
  case SC_SHMEM_WINDOW_PRESCAN:
#endif
    t8_shmem_allgather_common (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, intranode, internode);
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
