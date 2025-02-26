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

/** \file t8_forest_general.h
 * We define the forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_IO_H
#define T8_FOREST_IO_H

#include <t8_vtk.h>
T8_EXTERN_C_BEGIN ();

/* TODO: implement */
void
t8_forest_save (t8_forest_t forest);

/** Write the forest in a parallel vtu format. Extended version.
 * See \ref t8_forest_write_vtk for the standard version of this function.
 * Writes one master .pvtu file and each process writes in its own .vtu file.
 * If linked and not otherwise specified, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This may change in accordance with \a write_ghosts, \a write_curved and 
 * \a do_not_use_API, because the export of ghosts is not yet available with 
 * the VTK API and the export of curved elements is not available with the
 * inbuilt function to write ASCII files. The function will for example
 * still use the VTK API to satisfy \a write_curved, even if \a do_not_use_API 
 * is set to true.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 * \param [in]      forest              The forest to write.
 * \param [in]      fileprefix          The prefix of the files where the vtk will
 *                                      be stored. The master file is then fileprefix.pvtu
 *                                      and the process with rank r writes in the file
 *                                      fileprefix_r.vtu.
 * \param [in]      write_treeid        If true, the global tree id is written for each element.
 * \param [in]      write_mpirank       If true, the mpirank is written for each element.
 * \param [in]      write_level         If true, the refinement level is written for each element.
 * \param [in]      write_element_id    If true, the global element id is written for each element.
 * \param [in]      write_ghosts        If true, each process additionally writes its ghost elements.
 *                                      For ghost element the treeid is -1.
 * \param [in]      write_curved        If true, write the elements as curved element types from vtk.
 * \param [in]      do_not_use_API      Do not use the VTK API, even if linked and available.
 * \param [in]      num_data            Number of user defined double valued data fields to write.
 * \param [in]      data                Array of t8_vtk_data_field_t of length \a num_data
 *                                      providing the user defined per element data.
 *                                      If scalar and vector fields are used, all scalar fields
 *                                      must come first in the array.
 * \return  T8_SUBROUTINE_SUCCESS if successful, T8_SUBROUTINE_FAILED if not (process local).
 * See also \ref t8_forest_write_vtk .
 */
int
t8_forest_write_vtk_ext (t8_forest_t forest, const char *fileprefix, const int write_treeid, const int write_mpirank,
                         const int write_level, const int write_element_id, const int write_ghosts,
                         const int write_curved, int do_not_use_API, const int num_data, t8_vtk_data_field_t *data);

/** Write the forest in a parallel vtu format. Writes one master
 * .pvtu file and each process writes in its own .vtu file.
 * If linked, the VTK API is used.
 * If the VTK library is not linked, an ASCII file is written.
 * This function writes the forest elements, the tree id, element level, mpirank and element id as data.
 * Forest must be committed when calling this function.
 * This function is collective and must be called on each process.
 * For more options use \ref t8_forest_write_vtk_ext
 * \param [in]      forest              The forest to write.
 * \param [in]      fileprefix          The prefix of the files where the vtk will
 *                                      be stored. The master file is then fileprefix.pvtu
 *                                      and the process with rank r writes in the file
 *                                      fileprefix_r.vtu.
 * \return  T8_SUBROUTINE_SUCCESS if successful, T8_SUBROUTINE_FAILED if not (process local).
 */
int
t8_forest_write_vtk (t8_forest_t forest, const char *fileprefix);
T8_EXTERN_C_END ();

#endif /* !T8_FOREST_IO_H */
