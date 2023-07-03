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

#include "t8_vtk_parallel.hxx"

#if T8_WITH_VTK
#include <vtkXMLPDataReader.h>

vtk_read_success_t
t8_read_parallel (const char *filename, vtkSmartPointer < vtkDataSet > grid,
                  sc_MPI_Comm comm)
{
  /* Check if we can open the parallel file */
  FILE               *first_check;
  first_check = fopen (filename, "r");
  if (first_check == NULL) {
    t8_errorf ("Can not find the file %s\n", filename);
    return read_failure;
  }
  fclose (first_check);

  /* Setup parallel reader. */
  vtkSmartPointer < vtkXMLPDataReader > reader =
    vtkSmartPointer < vtkXMLPDataReader >::New ();
  if (!reader->CanReadFile (filename)) {
    return read_failure;
  }
  /* Get mpi size and rank */
  const int           total_num_pieces = reader->GetNumberOfPieces ();
  int                 mpiret;
  int                 mpisize;
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  int                 mpirank;
  mpiret = sc_MPI_Comm_rank (comm & mpirank);
  SC_CHECK_MPI (mpiret);

  /* Setup number of pieces to read on this proc. */
  int                 local_num_pieces = 0;
  if (mpisize >= total_num_pieces) {
    local_num_pieces = mpirank < mpisize ? 1 : 0;
  }
  else {
    local_num_pieces = total_num_pieces / mpisize * mpirank;
  }

  return read_failure;
}

#endif
