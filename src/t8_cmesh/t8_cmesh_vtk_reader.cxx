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

/** \file t8_cmesh_vtk_reader.cxx
* Implementation of a Reader for vtk/vtu files using the vtk-library.
* The functions can only be used when t8code is linked with the vtk-library.
*/

#include <t8_cmesh_vtk_reader.hxx>
#include "t8_cmesh_vtk_to_t8/t8_cmesh_vtk_unstructured.hxx"
#include "t8_cmesh_vtk_to_t8/t8_cmesh_vtk_polydata.hxx"

#if T8_WITH_VTK
#include <vtkDataObject.h>
#endif

T8_EXTERN_C_BEGIN ();

#if T8_WITH_VTK
int
t8_file_to_vtkGrid (const char *filename,
                    vtkDataSet * vtkGrid,
                    const int partition, const int main_proc,
                    sc_MPI_Comm comm, const vtk_file_type_t vtk_file_type)
{
  int                 main_proc_read_successful = 0;
  int                 mpirank;
  int                 mpiret;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (!partition || mpirank == main_proc) {
    switch (vtk_file_type) {
    case VTK_UNSTRUCTURED_FILE:
      t8_read_unstructured (filename, vtkGrid);
      break;
    case VTK_POLYDATA_FILE:
      t8_read_poly (filename, vtkGrid);
    default:
      vtkGrid = NULL;
      t8_errorf ("Filetype not supported.\n");
      break;
    }
    if (vtkGrid == NULL) {
      t8_errorf ("Could not read file\n");
      if (partition) {
        main_proc_read_successful = 0;
        sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc,
                      comm);
      }
    }
    main_proc_read_successful = 1;
  }
  if (partition) {
    sc_MPI_Bcast (&main_proc_read_successful, 1, sc_MPI_INT, main_proc, comm);
  }
  return main_proc_read_successful;
}
#endif

t8_cmesh_t
t8_cmesh_vtk_reader (const char *filename, const int partition,
                     const int main_proc, sc_MPI_Comm comm,
                     const vtk_file_type_t vtk_file_type)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh;
  int                 mpisize;
  int                 mpirank;
  int                 mpiret;
  /* Get the size of the communicator and the rank of the process. */
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Ensure that the main-proc is a valid proc. */
  T8_ASSERT (0 <= main_proc && main_proc < mpisize);
  T8_ASSERT (filename != NULL);
  int                 main_proc_read_successful = 0;

  vtkDataSet         *vtkGrid = NULL;
  T8_ASSERT (partition == 0 || (main_proc >= 0 && main_proc < mpisize));
  /* Read the file and set the pointer. */
  main_proc_read_successful =
    t8_file_to_vtkGrid (filename, vtkGrid, partition, main_proc, comm,
                        vtk_file_type);
  if (!main_proc_read_successful) {
    t8_global_errorf
      ("Main process (Rank %i) did not read the file successfully.\n",
       main_proc);
    t8_cmesh_destroy (&cmesh);
    return NULL;
  }
  else {
    /*TODO: translate the vtkGrid into a cmesh and partition it. */
    SC_ABORT
      ("The translation from a vtkGrid to a cmesh is not implemented yet.\n");
    return NULL;
  }
  T8_ASSERT (cmesh != NULL);
  return cmesh;
#else
  /* Return NULL if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
#endif
  return NULL;
}

T8_EXTERN_C_END ();
