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

#if T8_WITH_VTK
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#endif

T8_EXTERN_C_BEGIN ();
/**
 * Enumerator for all types of files readable by t8code. 
 */
typedef enum vtk_file_type
{
  VTK_FILE_ERROR = -1,          /*For Testing purpose. */
  VTK_UNSTRUCTURED_FILE = 0,
  VTK_POLYDATA_FILE = 1
} vtk_file_type_t;

/**
 * Translator between vtk-type of elements and t8code-elements. 
 * Not all elements are supported. Return T8_ECLASS_INVALID for unsupported
 * elements.  
 */
extern const t8_eclass_t t8_cmesh_vtk_type_to_t8_type[82] = {
  T8_ECLASS_INVALID, T8_ECLASS_VERTEX, T8_ECLASS_INVALID, T8_ECLASS_LINE,
  T8_ECLASS_INVALID, T8_ECLASS_TRIANGLE, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_TET,
  T8_ECLASS_HEX, T8_ECLASS_HEX, T8_ECLASS_PRISM, T8_ECLASS_PYRAMID,
  T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
  T8_ECLASS_INVALID, T8_ECLASS_INVALID
};

/**
 * Construct a cmesh given a filename refering to a vtk-file either containing an
 * unstructured grid, or vtk-polydata, most likely constructed by vtk. 
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
