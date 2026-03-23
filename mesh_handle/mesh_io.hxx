/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file mesh_io.hxx
 * In- and output of meshes. For example, writing a mesh to vtk format.
 */

#pragma once
#include <t8_forest/t8_forest_io.h>
namespace t8_mesh_handle
{

/** 
 * Write the mesh in a parallel vtu format. Extended version.
 * See \see write_mesh_to_vtk for the standard version of this function.
 * Writes one master .pvtu file and each process writes in its own .vtu file.
 * If linked and not otherwise specified, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This may change in accordance with \a write_ghosts, \a write_curved and 
 * \a do_not_use_API, because the export of ghosts is not yet available with 
 * the VTK API and the export of curved elements is not available with the
 * inbuilt function to write ASCII files. The function will for example
 * still use the VTK API to satisfy \a write_curved, even if \a do_not_use_API 
 * is set to true.
 * This function is collective and must be called on each process.
 * \param [in] mesh             The mesh to write.
 * \param [in] fileprefix       The prefix of the files where the vtk will be stored.
 *             The master file is then fileprefix.pvtu and the process with rank r writes in the file fileprefix_r.vtu.
 * \param [in] num_data         Number of user defined double valued data fields to write.
 * \param [in] data             Array of t8_vtk_data_field_t of length \a num_data providing the user defined 
 *              per element data. If scalar and vector fields are used, all scalar fields must come first in the array.
 * \param [in] write_treeid     If true, the global tree id of the underlying forest is written for each element.
 * \param [in] write_mpirank    If true, the mpirank is written for each element.
 * \param [in] write_level      If true, the refinement level is written for each element.
 * \param [in] write_element_id If true, the global element id is written for each element.
 * \param [in] write_ghosts     If true, each process additionally writes its ghost elements.
 *              For ghost element the treeid of the underlying forest is -1.
 * \param [in] write_curved     If true, write the elements as curved element types from vtk.
 * \param [in] do_not_use_API   Do not use the VTK API, even if linked and available.
 * \return  True if successful, false if not (process local).
 */
template <typename TMeshClass>
int
write_mesh_to_vtk_ext (TMeshClass &mesh, const char *fileprefix, const int num_data, t8_vtk_data_field_t *data,
                       bool write_treeid = false, bool write_mpirank = true, bool write_level = true,
                       bool write_element_id = true, bool write_ghosts = false, bool write_curved = false,
                       bool do_not_use_API = false)
{
  return t8_forest_write_vtk_ext (mesh->get_forest (), fileprefix, write_treeid, write_mpirank, write_level,
                                  write_element_id, write_ghosts, write_curved, do_not_use_API, num_data, data);
}

/** 
 * Write the mesh in a parallel vtu format. Writes one master
 * .pvtu file and each process writes in its own .vtu file.
 * If linked, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This function writes the elements, level, mpirank, element id and the tree id of the underlying forest as data.
 * This function is collective and must be called on each process.
 * For more options use \see write_mesh_to_vtk_ext.
 * \param [in] mesh       The mesh to write.
 * \param [in] fileprefix The prefix of the files where the vtk will be stored.
 *             The master file is then fileprefix.pvtu and the process with rank r writes in the file fileprefix_r.vtu.
 * \return  True if successful, false if not (process local).
 */
template <typename TMeshClass>
int
write_mesh_to_vtk (TMeshClass &mesh, const char *fileprefix)
{
  return t8_forest_write_vtk (mesh->get_forest (), fileprefix);
}

}  // namespace t8_mesh_handle
