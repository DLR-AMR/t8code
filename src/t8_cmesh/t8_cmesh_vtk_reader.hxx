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

#ifndef T8_CMESH_VTK_READER
#define T8_CMESH_VTK_READER

#include <t8_cmesh.h>

T8_EXTERN_C_BEGIN ();
#if T8_WITH_VTK
/* look_up Table to transform vtkCellType into T8_ECLASS
   T8_ECLASS_ZERO, if a the vtkCellType is not supported by t8code.
   see https://vtk.org/doc/nightly/html/vtkCellType_8h.html to check.*/
const t8_eclass_t   t8_cmesh_vtk_type_to_t8_type[82] = {
  T8_ECLASS_ZERO, T8_ECLASS_VERTEX, T8_ECLASS_ZERO, T8_ECLASS_LINE,
  T8_ECLASS_ZERO, T8_ECLASS_TRIANGLE, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_QUAD, T8_ECLASS_TET, T8_ECLASS_ZERO,
  T8_ECLASS_HEX, T8_ECLASS_PRISM, T8_ECLASS_PYRAMID, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO, T8_ECLASS_ZERO,
  T8_ECLASS_ZERO, T8_ECLASS_ZERO
};
#endif

/**
 * Construct a cmesh given a filename and a number of files (for parallel reader)
 * CAREFULL: This is in production and this header will probably change! Update
 * as the function progresses
 * 
 * \param[in] filename      The name of the file 
 * \param[in] num_files     Number of files to read from
 * @return t8_cmesh_t       The cmesh described by the files
 */
t8_cmesh_t          t8_cmesh_read_from_vtk (const char *filename,
                                            const int num_files,
                                            sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* T8_CMESH_VTK_READER */
