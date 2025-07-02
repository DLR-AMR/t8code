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

#ifndef T8_VTK_READER
#define T8_VTK_READER

#include <t8_cmesh.h>
#include <t8_vtk/t8_vtk_types.h>

#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>

/**
 * Given a pointer to a vtkDataSet a cmesh representing the vtkDataSet is
 * constructed and can be shared over the processes. The vtk data arrays will be associated with the trees.
 * 
 * \param[in] vtkGrid A pointer to a vtkDataSet
 * \param[in] partition Flag if the cmesh should be partitioned
 * \param[in] main_proc The main reading process
 * \param[in] distributed_grid Flag if the vtkGrid is distributed over several procs. 
 * \param[in] comm The communicator. 
 * \param[in] package_id The package id of the application. It is generated with the usage of \ref sc_package_register.
 * \param[in] starting_key If the application already registered attributes, the starting key is used so that the existing attributes are not overwritten.
 * \return The committed cmesh 
 */
t8_cmesh_t
t8_vtkGrid_to_cmesh (vtkSmartPointer<vtkDataSet> vtkGrid, const int partition, const int main_proc,
                     const int distributed_grid, sc_MPI_Comm comm, const int package_id, const int starting_key);

/**
 * Given a pointer to a vtkDataSet a vtkPointSet storing a set of points of
 * is constructed. The cell data of vtkDataSet is mapt on the points of vtkPointSet. 
 * 
 * \param[in] vtkGrid A pointer to a vtkDataSet
 * \return A pointer to a vtkPointSet 
 */
vtkSmartPointer<vtkPointSet>
t8_vtkGrid_to_vtkPointSet (vtkSmartPointer<vtkDataSet> vtkGrid);

/**
 * Given a filename to a vtkUnstructuredGrid or vtkPolyData read the file and
 * construct a abstract class to specify dataset behavior. The file is read and
 * stored in a vtkDataSet.
 * \note This function is only available if t8code is linked against VTK. 
 * 
 * \param[in] filename      The name of the file
 * \param[in] partition     Flag if the constructed mesh should be partitioned
 * \param[in] main_proc     The main reading processor
 * \param[in] comm          An mpi-communicator
 * \param[in] vtk_file_type A vtk-filetype that is readable by t8code. 
 * \return                  Pointer to vtkDataSet
 */
vtkSmartPointer<vtkDataSet>
t8_vtk_reader (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
               const vtk_file_type_t vtk_file_type);

/**
 * Given a filename to a vtkUnstructuredGrid or vtkPolyData read the file and
 * a set of points is constructed. This is a two stage process. First the file
 * is read and stored in a vtkDataSet using \a t8_vtk_reader and 
 * \a t8_file_to_vtkGrid. In the second stage a vtkPointSet is constructed from 
 * the vtkDataSet using \a t8_vtkGrid_to_vtkPointSet. 
 * 
 * Both stages use the vtk-library, therefore the function is only available if 
 * t8code is linked against VTK. 
 * 
 * \param[in] filename      The name of the file
 * \param[in] partition     Flag if the constructed mesh should be partitioned
 * \param[in] main_proc     The main reading processor
 * \param[in] comm          An mpi-communicator
 * \param[in] vtk_file_type A vtk-filetype that is readable by t8code. 
 * \return                  Pointer to vtkDataSet      
 */
vtkSmartPointer<vtkPointSet>
t8_vtk_reader_pointSet (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
                        const vtk_file_type_t vtk_file_type);
/**
 * Given a filename to a vtkUnstructuredGrid or vtkPolyData read the file and
 * construct a cmesh. This is a two stage process. First the file is read and
 * stored in a vtkDataSet using \a t8_vtk_reader and \a t8_file_to_vtkGrid. 
 * In the second stage a cmesh is constructed from the vtkDataSet using \a t8_vtkGrid_to_cmesh.
 * The vtk data arrays will be associated with the cmesh trees and saved as tree attributes using the 
 * provided user application \a package_id. The keys for the tree attributes start at \a starting_key.
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
t8_vtk_reader_cmesh (const char *filename, const int partition, const int main_proc, sc_MPI_Comm comm,
                     const vtk_file_type_t vtk_file_type, const int package_id, const int starting_key);

#endif /* T8_VTK_READER */
