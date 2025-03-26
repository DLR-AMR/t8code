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
#ifndef T8_WITH_VTK

#include <t8_vtk/t8_vtk_reader.hxx>
#include <t8_vtk/t8_vtk_types.h>
#include <t8_cmesh.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>

T8_EXTERN_C_BEGIN ();

/**
 * If the vertices of a tree describe a negative \param, 
 * permute the tree vertices. 
 * 
 * \param[in, out] tree_vertices The vertices of a tree
 * \param[in] eclass             The eclass of the tree.
 */

void
t8_vtk_reader_pointSet ([[maybe_unused]] const char *filename, [[maybe_unused]] const int partition,
                        [[maybe_unused]] const int main_proc, [[maybe_unused]] sc_MPI_Comm comm,
                        [[maybe_unused]] const vtk_file_type_t vtk_file_type)
{
  /* Return NULL since not linked against vtk */
  t8_global_errorf (
    "WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
}

t8_cmesh_t
t8_vtk_reader_cmesh ([[maybe_unused]] const char *filename, [[maybe_unused]] const int partition,
                     [[maybe_unused]] const int main_proc, [[maybe_unused]] sc_MPI_Comm comm,
                     [[maybe_unused]] const vtk_file_type_t vtk_file_type)
{
  /* Return NULL since not linked against vtk */
  t8_global_errorf (
    "WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return NULL;
}

T8_EXTERN_C_END ();

#endif
