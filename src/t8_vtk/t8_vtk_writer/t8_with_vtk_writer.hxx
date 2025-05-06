/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#include <t8_cmesh.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_vtk/t8_vtk_writer_helper.hxx>
#include <t8_vtk/t8_vtk_write_ASCII.hxx>

#include <string>
#include <t8_vtk.h>
#include <t8_types/t8_vec.hxx>

template<typename grid_t>
class vtk_enabled_writer : public vtk_writer<grid_t>{

 public:
 void t8_grid_tree_to_vtk_cells(const t8_forest_t forest, [[maybe_unused]] vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                                vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
                                vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
                                const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
                                const bool ghosts, const t8_locidx_t itree)

void t8_grid_tree_to_vtk_cells (const t8_cmesh_t cmesh, [[maybe_unused]] vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
                                vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
                                vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
                                [[maybe_unused]] const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id,
                                const t8_gloidx_t offset, const bool ghosts, const t8_locidx_t itree)
}