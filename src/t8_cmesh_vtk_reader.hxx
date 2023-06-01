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

/** \file t8_cmesh_vtk_reader.hxx
* Header for the vtk-reader. 
*/

#ifndef T8_CMESH_VTK_READER
#define T8_CMESH_VTK_READER

#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_vtk_to_t8/t8_vtk_types.h"

#if T8_WITH_VTK
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#endif

T8_EXTERN_C_BEGIN ();

#if T8_WITH_VTK
/**
 * Read a vtk-file and ShallowCopy its content into a vtkDataSet.
 * The success (or failure) of the reading process is communicated 
 * over all processes.
 * 
 * \param[in] filename      The name of the file to read
 * \param[in, out] vtkGrid  A pointer to a vtkDataSet. We ShallowCopy the grid there.
 * \param[in] partition     Flag if the input is read partitioned
 * \param[in] main_proc     The main reading proc.
 * \param[in] comm          A communicator.
 * \param[in] vtk_file_type The type of the Data in the file.
 * \return                  0 if the file was read successfully, 1 otherwise.                
 */
vtk_read_success_t  t8_file_to_vtkGrid (const char *filename,
                                        vtkSmartPointer < vtkDataSet >
                                        vtkGrid, const int partition,
                                        const int main_proc, sc_MPI_Comm comm,
                                        const vtk_file_type_t vtk_file_type);

/**
 * Given a pointer to a vtkDataSet a cmesh representing the vtkDataSet is
 * constructed and can be shared over the processes. 
 * 
 * \param[in] vtkGrid A pointer to a vtkDataSet
 * \param[in] partition Flag if the cmesh should be partitioned
 * \param[in] main_proc The main reading process
 * \param[in] comm The communicator. 
 * \return t8_cmesh_t 
 */
t8_cmesh_t          t8_vtkGrid_to_cmesh (vtkSmartPointer < vtkDataSet >
                                         vtkGrid, const int partition,
                                         const int main_proc,
                                         sc_MPI_Comm comm);

#endif

/**
 * Given a filename to a vtkUnstructuredGrid or vtkPolyData read the file and
 * construct a cmesh. This is a two stage process. First the file is read and
 * stored in a vtkDataSet using \a t8_file_to_vtkGrid. 
 * In the second stage a cmesh is constructed from the vtkDataSet using \a t8_vtkGrid_to_cmesh. 
 * 
 * Both stages use the vtk-library, therefore the function is only available if 
 * t8code is linked against VTK. 
 * 
 * 
 * \param[in] filename      The name of the file
 * \param[in] partition     Flag if the constructed mesh should be partitioned
 * \param[in] main_proc     The main reading processor
 * \param[in] comm          An mpi-communicator
 * \param[in] vtk_file_type A vtk-filetype that is readable by t8code. 
 * \return                  A commited cmesh.       
 */
t8_cmesh_t          t8_cmesh_vtk_reader (const char *filename,
                                         const int partition,
                                         const int main_proc,
                                         sc_MPI_Comm comm,
                                         const vtk_file_type_t vtk_file_type);

T8_EXTERN_C_END ();

#endif /* T8_CMESH_VTK_READER */
