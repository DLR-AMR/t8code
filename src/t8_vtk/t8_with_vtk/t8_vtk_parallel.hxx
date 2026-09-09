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

/**
 * \file t8_vtk_parallel.hxx
 * This file contains all helper functions to translate parallel vtk files. 
 */

#ifndef T8_VTK_PARALLEL_HXX
#define T8_VTK_PARALLEL_HXX

#include <t8.h>
#include <t8_vtk/t8_vtk_types.h>

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

/**
 * Given a filename to a parallel vtk file (for example .pvtu) and its data files, 
 * read a piece of the data files (like .vtu, ...).
 * 
 * \param[in] filename  The name of a parallel vtk file to a distributed vtkUnstructuredGrid
 * \param[in] grid On input a vtkSmartPointer, that will hold the grid described
 *                 by the pieces read on this proc.
 * \param[in] comm The MPI communicator to use for reading the file.
 * \returns        non-zero on success, zero if the reading failed. 
 */
vtk_read_success_t
t8_read_parallel_unstructured (const char *filename, vtkSmartPointer<vtkDataSet> grid, sc_MPI_Comm comm);

/**
 * Given a filename to a parallel vtk file (for example .pvtp) and its data files, 
 * read a piece of the data files (like  .vtp, ...).
 * 
 * \param[in] filename  The name of a parallel vtk file to a distributed vtkPolyData
 * \param[in] grid On input a vtkSmartPointer, that will hold the grid described
 *                 by the pieces read on this proc.
 * \param[in] comm The MPI communicator to use for reading the file.
 * \returns        read_success on success, read_failure if the reading failed.
 */
vtk_read_success_t
t8_read_parallel_polyData (const char *filename, vtkSmartPointer<vtkDataSet> grid, sc_MPI_Comm comm);

#endif /* T8_VTK_PARALLEL_HXX */
