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

#ifndef T8_VTK_WRITER_HXX
#define T8_VTK_WRITER_HXX

#include <t8_cmesh.h>
#include <t8_vtk.h>
#include <t8.h>
#include "t8_vtk/t8_vtk_writer_impl.hxx"

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#endif

/*
 TODO: make this a class so we dont have to trace what we want to write out

template <typename grid_t>
class vtk_writer
{
public:

private:

    bool write_treeid;
    bool write_mpirank;
    bool write_level;
    bool write_element_id;
    bool write_ghosts;
    bool curved_flag;
    int num_data;
    t8_vtk_data_field_t *data;
    sc_MPI_Comm comm;
}; */

/**
 * \brief 
 * 
 * \param[in] grid 
 * \param[in] fileprefix 
 * \param[in] write_treeid 
 * \param[in] write_mpirank 
 * \param[in] write_level 
 * \param[in] write_element_id 
 * \param[in] curved_flag 
 * \param[in] write_ghosts 
 * \param[in] num_data 
 * \param[in] data 
 * \return true 
 * \return false 
 */
template <typename grid_t>
int
t8_write_vtk_via_API (const grid_t grid, std::string fileprefix, const int write_treeid, const int write_mpirank,
                      const int write_level, const int write_element_id, const int curved_flag, const int write_ghosts,
                      const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
  return t8_write_vtk (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id, curved_flag,
                       write_ghosts, num_data, data, comm);
}

#endif /* T8_VTK_WRITER_HXX */