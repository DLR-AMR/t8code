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

/**
 * This file contains all helper functions to translate vtk-polydata. 
 * 
 */

#ifndef T8_CMESH_VTK_POLYDATA
#define T8_CMESH_VTK_POLYDATA

#include <t8.h>
#include "t8_vtk_types.h"
#if T8_ENABLE_VTK
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>

/**
 * Given a filename to a file containing Polydata, read
 * the file using the vtk-library. 
 * 
 * \param[in] filename  The name of the file
 * \param[in, out] grid On input a vtkSmartPointer, that will hold the grid described in
 *                      \a filename.
 * \returns             non-zero on success, zero if the reading failed.              
 * 
 */
vtk_read_success_t
t8_read_polyData (const char *filename, vtkDataSet *grid);

#endif /* T8_ENABLE_VTK */
#endif /* T8_CMESH_VTK_POLYDATA */
