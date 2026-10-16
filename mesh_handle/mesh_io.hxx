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
#include <t8_vector/t8_vec.hxx>
#include "concepts.hxx"
#include <span>
#include <string>
#include <type_traits>
#include <vector>

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
 * built-in function to write ASCII files. The function will for example
 * still use the VTK API to satisfy \a write_curved, even if \a do_not_use_API 
 * is set to true.
 * This function is collective and must be called on each process.
 * \param [in] mesh             The mesh to write.
 * \param [in] fileprefix       The prefix of the files where the vtk will be stored.
 *             The master file is then fileprefix.pvtu and the process with rank r writes in the file fileprefix_r.vtu.
 * \param [in] num_data         Number of user-defined double valued data fields to write.
 * \param [in] data             Array of t8_vtk_data_field_t of length \a num_data providing the user-defined 
 *              per element data. If scalar and vector fields are used, all scalar fields must come first in the array.
 * \param [in] write_treeid     If true, the global tree id of the underlying forest is written for each element.
 * \param [in] write_mpirank    If true, the mpirank is written for each element.
 * \param [in] write_level      If true, the refinement level is written for each element.
 * \param [in] write_element_id If true, the global element id is written for each element.
 * \param [in] write_ghosts     If true, each process additionally writes its ghost elements.
 *              For ghost element the treeid of the underlying forest is -1.
 * \param [in] write_curved     If true, write the elements as curved element types from vtk.
 * \param [in] do_not_use_API   If true, do not use the VTK API, even if linked and available.
 * \return  True if successful, false if not (process local).
 */
template <T8MeshType TMeshClass>
int
write_mesh_to_vtk_ext (TMeshClass &mesh, const char *fileprefix, const int num_data, t8_vtk_data_field_t *data,
                       bool write_treeid = false, bool write_mpirank = true, bool write_level = true,
                       bool write_element_id = true, bool write_ghosts = false, bool write_curved = false,
                       bool do_not_use_API = false)
{
  return t8_forest_write_vtk_ext (mesh.get_forest (), fileprefix, write_treeid, write_mpirank, write_level,
                                  write_element_id, write_ghosts, write_curved, do_not_use_API, num_data, data);
}

/** A named per element data field for the vtk output, see \ref write_mesh_to_vtk_ext_data.
 *  \tparam TValue The value type per element: double for scalar fields and t8_3D_vec for vector fields.
 */
template <typename TValueType>
struct vtk_data_field
{
  std::string name;                   /**< The name of the field as it should be written to the file. */
  std::span<const TValueType> values; /**< One value per local element, optionally also one per ghost element. */
};
using vtk_scalar_field = vtk_data_field<double>;    /**< One double value per element (T8_VTK_SCALAR). */
using vtk_vector_field = vtk_data_field<t8_3D_vec>; /**< One 3D vector value per element (T8_VTK_VECTOR). */

namespace detail
{
/** Helper function to convert a \ref vtk_data_field into the t8_vtk_data_field_t expected by the vtk writer.
 *  The values are copied into \a buffer, since the vtk writer expects a flat double array.
 *  This is used in \ref write_mesh_to_vtk_ext_data. 
 *  \tparam TMeshClass The mesh handle class.
 *  \tparam TValueType The type of the values in the vtk data field.
 *  \param [in] mesh The mesh the data field belongs to.
 *  \param [in] field The vtk data field to convert.
 *  \param [in] write_ghosts If true, ghost elements are written too.
 *  \param [out] buffer The buffer for the values. It has to stay alive until the vtk output is written.
 *  \param [out] vtk_field The converted data field pointing to \a buffer. 
*/
template <T8MeshType TMeshClass, typename TValueType>
void
to_vtk_data_field (const TMeshClass &mesh, const vtk_data_field<TValueType> &field, bool write_ghosts,
                   std::vector<double> &buffer, t8_vtk_data_field_t &vtk_field)
{
  constexpr bool is_scalar = std::is_same_v<TValueType, double>;
  const t8_locidx_t num_local_elements = mesh.get_num_local_elements ();
  const t8_locidx_t num_ghosts = mesh.get_num_ghosts ();
  T8_ASSERT (field.values.size () == static_cast<size_t> (num_local_elements)
             || field.values.size () == static_cast<size_t> (num_local_elements + num_ghosts));
  for (const TValueType &value : field.values) {
    if constexpr (is_scalar) {
      buffer.push_back (value);
    }
    else {
      buffer.insert (buffer.end (), value.begin (), value.end ());
    }
  }
  /** The vtk writer expects exactly one value per cell. If ghosts are written, but no values are given, 
    * all of them are set to zero. If ghosts are not written, the values are dropped. 
  */
  const t8_locidx_t num_cells = write_ghosts ? num_local_elements + num_ghosts : num_local_elements;
  buffer.resize (static_cast<size_t> (num_cells) * (is_scalar ? 1 : 3), 0.0);

  vtk_field.type = is_scalar ? T8_VTK_SCALAR : T8_VTK_VECTOR;
  snprintf (vtk_field.description, BUFSIZ, "%s", field.name.c_str ());
  vtk_field.data = buffer.data ();
}
} /* Namespace detail */

/**
 * Write the mesh in a parallel vtu format together with user-defined per element data fields.
 * Convenience wrapper around \ref write_mesh_to_vtk_ext that takes the data as named std::vectors,
 * so there is no need to allocate arrays or to fill t8_vtk_data_field_t structures manually. 
 * \param [in] mesh       The mesh to write.
 * \param [in] fileprefix The prefix of the files where the vtk will be stored.
 *             The master file is then fileprefix.pvtu and the process with rank r writes in the file fileprefix_r.vtu.
 * \param [in] scalar_fields  The named scalar per element data fields to write.
 * \param [in] vector_fields  The named vector per element data fields to write.
 * \param [in] write_treeid  Whether to write the tree id of the underlying forest.
 * \param [in] write_mpirank Whether to write the MPI rank of the process owning the element.
 * \param [in] write_level   Whether to write the level of the element.
 * \param [in] write_element_id Whether to write the element id.
 * \param [in] write_ghosts  Whether to write ghost elements.
 * \param [in] write_curved  Whether to write curved elements.
 * \param [in] do_not_use_API Whether to bypass the VTK API and write an ASCII file directly.
*/
template <T8MeshType TMeshClass>
int
write_mesh_to_vtk_ext_data (TMeshClass &mesh, const char *fileprefix,
                            const std::vector<vtk_scalar_field> &scalar_fields,
                            const std::vector<vtk_vector_field> &vector_fields = {}, bool write_treeid = false,
                            bool write_mpirank = true, bool write_level = true, bool write_element_id = true,
                            bool write_ghosts = false, bool write_curved = false, bool do_not_use_API = false)
{
  const int num_data = static_cast<int> (scalar_fields.size () + vector_fields.size ());
  std::vector<std::vector<double>> buffers (num_data);
  std::vector<t8_vtk_data_field_t> vtk_data (num_data);
  /* The vtk writer expects all scalar fields before vector fields. */
  int idata = 0;
  for (const vtk_scalar_field &field : scalar_fields) {
    detail::to_vtk_data_field (mesh, field, write_ghosts, buffers[idata], vtk_data[idata]);
    ++idata;
  }
  for (const vtk_vector_field &field : vector_fields) {
    detail::to_vtk_data_field (mesh, field, write_ghosts, buffers[idata], vtk_data[idata]);
    ++idata;
  }
  return write_mesh_to_vtk_ext (mesh, fileprefix, num_data, vtk_data.data (), write_treeid, write_mpirank, write_level,
                                write_element_id, write_ghosts, write_curved, do_not_use_API);
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
template <T8MeshType TMeshClass>
int
write_mesh_to_vtk (TMeshClass &mesh, const char *fileprefix)
{
  return t8_forest_write_vtk (mesh.get_forest (), fileprefix);
}

}  // namespace t8_mesh_handle
