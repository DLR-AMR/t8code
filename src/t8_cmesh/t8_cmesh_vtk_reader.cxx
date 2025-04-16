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

T8_EXTERN_C_BEGIN ();

t8_cmesh_t
t8_cmesh_vtk_reader (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
                     const vtk_file_type_t vtk_file_type, const int package_id, const int starting_key)
{
  return t8_vtk_reader_cmesh (filename, partition, main_proc, comm, vtk_file_type, package_id, starting_key);
}

T8_EXTERN_C_END ();
