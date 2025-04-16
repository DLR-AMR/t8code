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

#include <t8_vtk/t8_vtk_reader.hxx>

T8_EXTERN_C_BEGIN ();

/**
 * Given a filename to a vtkUnstructuredGrid or vtkPolyData read the file and
 * construct a cmesh. This is a two stage process. First the file is read and
 * stored in a vtkDataSet using \a t8_vtk_reader and \a t8_file_to_vtkGrid. 
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
 * \param[in] package_id    The package id of the application. It is generated with the usage of \ref sc_package_register.  
 * \param[in] starting_key  If the application already registered attributes, the starting key is used so that the existing attributes are not overwritten.
 * \return                  A committed cmesh.
 */
t8_cmesh_t
t8_cmesh_vtk_reader (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
                     const vtk_file_type_t vtk_file_type, const int package_id, const int starting_key);

T8_EXTERN_C_END ();

#endif /* T8_CMESH_VTK_READER */
