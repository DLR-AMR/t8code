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

#endif
