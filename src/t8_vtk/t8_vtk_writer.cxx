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

#include <t8_vtk/t8_vtk_writer.hxx>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#endif

#if T8_WITH_VTK
/**
 * \brief template specialization for forests. 
 * 
 */
template <>
void
vtk_writer<t8_forest_t>::t8_grid_tree_to_vtk_cells (
  const t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
  const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
  const bool ghosts, const t8_locidx_t itree)
{
  /* For both ghosts and pure-local trees iterate over all elements and translate them into a vtk cell. */
  if (ghosts) {
    const t8_locidx_t num_ghosts = t8_forest_ghost_tree_num_elements (forest, itree);
    for (t8_locidx_t ielem_ghost = 0; ielem_ghost < num_ghosts; ielem_ghost++) {
      const t8_element_t *element = t8_forest_ghost_get_element (forest, itree, ielem_ghost);
      this->t8_grid_element_to_vtk_cell (forest, element, itree + num_local_trees, offset, true, *elem_id, point_id,
                                         cellTypes, points, cellArray, vtk_treeid, vtk_mpirank, vtk_level,
                                         vtk_element_id);
      (*elem_id)++;
    }
  }
  else {
    const t8_locidx_t elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    /* We iterate over all elements in the tree */
    for (t8_locidx_t ielement = 0; ielement < elems_in_tree; ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      this->t8_grid_element_to_vtk_cell (forest, element, itree, offset, true, *elem_id, point_id, cellTypes, points,
                                         cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);
      (*elem_id)++;

    } /* end of loop over elements */
  }
  return;
}

/**
 * \brief template specialization for cmeshes. 
 * 
 */
template <>
void
vtk_writer<t8_cmesh_t>::t8_grid_tree_to_vtk_cells (
  const t8_cmesh_t cmesh, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_treeid, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_mpirank,
  vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_level, vtkSmartPointer<t8_vtk_gloidx_array_type_t> vtk_element_id,
  vtkSmartPointer<vtkCellArray> cellArray, vtkSmartPointer<vtkPoints> points, int *cellTypes,
  const t8_locidx_t num_local_trees, t8_gloidx_t *elem_id, long int *point_id, const t8_gloidx_t offset,
  const bool ghosts, const t8_locidx_t itree)
{
  /* A cmesh does not have any further elements, we can call the translatore directly. */
  this->t8_grid_element_to_vtk_cell (cmesh, NULL, itree, offset, ghosts, *elem_id, point_id, cellTypes, points,
                                     cellArray, vtk_treeid, vtk_mpirank, vtk_level, vtk_element_id);
  (*elem_id)++;
  return;
}
#endif /* T8_WITH_VTK */

template <>
bool
vtk_writer<t8_forest_t>::write_ASCII (const t8_forest_t forest)
{
  return t8_forest_vtk_write_ASCII (forest, this->fileprefix.c_str (), this->write_treeid, this->write_mpirank,
                                    this->write_level, this->write_element_id, this->write_ghosts, this->num_data,
                                    this->data);
}

template <>
bool
vtk_writer<t8_cmesh_t>::write_ASCII (const t8_cmesh_t forest)
{
  return t8_cmesh_vtk_write_ASCII (forest, this->fileprefix.c_str ());
}

/* Implementation of the c-interface */
T8_EXTERN_C_BEGIN ();

int
t8_forest_vtk_write_file_via_API (const t8_forest_t forest, const char *fileprefix, const int write_treeid,
                                  const int write_mpirank, const int write_level, const int write_element_id,
                                  const int curved_flag, const int write_ghosts, const int num_data,
                                  t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_with_API (forest);
}

int
t8_forest_vtk_write_file (const t8_forest_t forest, const char *fileprefix, const int write_treeid,
                          const int write_mpirank, const int write_level, const int write_element_id, int write_ghosts,
                          const int num_data, t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, false,
                                  std::string (fileprefix), num_data, data, t8_forest_get_mpicomm (forest));
  return writer.write_ASCII (forest);
}

int
t8_cmesh_vtk_write_file_via_API (const t8_cmesh_t cmesh, const char *fileprefix, sc_MPI_Comm comm)
{
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), comm);
  return writer.write_with_API (cmesh);
}

int
t8_cmesh_vtk_write_file (const t8_cmesh_t cmesh, const char *fileprefix)
{
  /* No mpi Communicator is needed for ASCII output*/
  vtk_writer<t8_cmesh_t> writer (std::string (fileprefix), sc_MPI_COMM_NULL);
  return writer.write_ASCII (cmesh);
}

#if T8_WITH_VTK
void
t8_forest_to_vtkUnstructuredGrid (const t8_forest_t forest, vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid,
                                  const int write_treeid, const int write_mpirank, const int write_level,
                                  const int write_element_id, const int write_ghosts, const int curved_flag,
                                  const int num_data, t8_vtk_data_field_t *data)
{
  vtk_writer<t8_forest_t> writer (write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, curved_flag,
                                  std::string (""), num_data, data, t8_forest_get_mpicomm (forest));
  writer.grid_to_vtkUnstructuredGrid (forest, unstructuredGrid);
}
#endif

T8_EXTERN_C_END ();
