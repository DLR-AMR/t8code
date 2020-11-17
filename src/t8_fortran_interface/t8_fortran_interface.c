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

#include <t8_cmesh.h>
#include <t8_fortran_interface/t8_fortran_interface.h>

/* Wrapper around sc_init (...), t8_init (...) */
void
t8_fortran_init_all_ (sc_MPI_Comm * comm)
{
#if T8_ENABLE_MPI
  /* Initialize sc */
  printf ("Sc init Start\n");
  sc_init (*comm, 1, 1, NULL, SC_LP_DEFAULT);
  printf ("Sc init\n");
  /* Initialize t8code */
  t8_init (SC_LP_DEFAULT);
  printf ("t8 init\n");
#else
  SC_ABORT ("t8code was not configured with MPI support.")
#endif
}

void
t8_fortran_init_all (sc_MPI_Comm * comm)
{
  sc_MPI_Comm         comm2 = sc_MPI_COMM_WORLD;
  int                 rank;

  printf ("Init all with comm %lu\n", (long unsigned) comm);
  sc_MPI_Comm_rank (*comm, &rank);
  printf ("rank = %i\n", rank);

  t8_fortran_init_all_ (&comm2);
}

/* Wrapper around sc_finalize */
void
t8_fortran_finalize ()
{
  sc_finalize ();
}

/* Build C MPI comm from Fortran MPI Comm. */
sc_MPI_Comm        *
t8_fortran_MPI_Comm_new (MPI_Fint Fcomm)
{
#if T8_ENABLE_MPI
  sc_MPI_Comm        *Ccomm = malloc (sizeof (*Ccomm));
  *Ccomm = MPI_Comm_f2c (Fcomm);
  printf ("Created comm %lu\n", (long unsigned) Ccomm);
  return Ccomm;
#else
  SC_ABORT ("t8code was not configured with MPI support.")
    return NULL;
#endif
}

/* Delete C MPI Comm. */
void
t8_fortran_MPI_Comm_delete (sc_MPI_Comm * Ccomm)
{
#if T8_ENABLE_MPI
  free (Ccomm);
#else
  SC_ABORT ("t8code was not configured with MPI support.")
#endif
}

t8_cmesh_t
t8_cmesh_new_periodic_tri_wrap (sc_MPI_Comm * Ccomm)
{
  return t8_cmesh_new_periodic_tri (*Ccomm);
}
