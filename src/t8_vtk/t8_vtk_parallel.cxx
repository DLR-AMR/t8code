/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPDataReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkAppendFilter.h>

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

  char                tmp[BUFSIZ], *extension;
  strcpy (tmp, filename);
  extension = strrchr (tmp, '.') + 1;
  T8_ASSERT (strcmp (extension, ""));

  /* Setup parallel reader. */
  vtkSmartPointer < vtkXMLPDataReader > reader = NULL;

  if (strcmp (extension, "pvtu") == 0) {
    reader = vtkSmartPointer < vtkXMLPUnstructuredGridReader >::New ();
  }
  else if (strcmp (extension, "pvtu") == 0) {
    reader = vtkSmartPointer < vtkXMLPPolyDataReader >::New ();
  }

  if (!reader->CanReadFile (filename)) {
    t8_errorf ("Unable to read file.\n");
    return read_failure;
  }

  reader->SetFileName (filename);
  reader->UpdateInformation ();

  /*  Get the number of files to read. */
  const int           total_num_pieces = reader->GetNumberOfPieces ();
  /* Get mpi size and rank */
  int                 mpiret;
  int                 mpisize;
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  int                 mpirank;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Setup number of pieces to read on this proc. */
  int                 last_piece = -1;
  int                 first_piece = 0;
  if (mpisize >= total_num_pieces) {
    /* The first n-procs read a piece each. */
    first_piece = mpirank;
    last_piece = first_piece + ((mpirank < mpisize) ? 1 : 0);
  }
  else {
    first_piece = total_num_pieces / mpisize * mpirank;
    last_piece = first_piece;
    const int           prev_proc_first_piece =
      total_num_pieces / mpisize * (mpirank - 1);
    last_piece +=
      first_piece == prev_proc_first_piece ? 0 : total_num_pieces / mpisize;
    if (first_piece == total_num_pieces / mpisize * (mpisize - 1)) {
      /* Read the last chunk of data */
      last_piece = total_num_pieces - first_piece;
    }
  }

  /* Read the pieces if there are any pieces to read on this proc. */
  if (first_piece < last_piece) {
    vtkNew < vtkAppendFilter > append;
    for (int ipiece = first_piece; ipiece < last_piece; ipiece++) {
      reader->UpdatePiece (ipiece, total_num_pieces, 0);
      append->AddInputData (reader->GetOutput ());
    }
    /* Merge all read grids together */
    append->Update ();
    append->MergePointsOn ();
    grid->ShallowCopy (append->GetOutputAsDataSet ());
  }
  else {
    /* Initialize the grid, but don't construct any cells. 
     * simplifies further processing of the grid on multiple procs. */
    grid->Initialize ();
  }
  return read_success;
}

#endif
