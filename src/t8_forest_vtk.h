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

/** file t8_forest_vtk.h
 */

/* TODO: Document this file */

#ifndef T8_FOREST_VTK_H
#define T8_FOREST_VTK_H

#include <t8_vtk.h>
#include <t8_forest.h>
/*
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkPyramid.h>
#include <vtkWedge.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkUnsignedCharArray.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_vtk.h>
#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"*/

T8_EXTERN_C_BEGIN ();
/* function declarations */

void
 
             t8_write_vtk_via_API (t8_forest_t forest, const char *filename);

/** Write the forest in .pvtu file format. Writes one .vtu file per
 * process and a meta .pvtu file.
 * \param [in]  forest    The forest.
 * \param [in]  fileprefix  The prefix of the output files.
 * \param [in]  write_treeid If true, the global tree id is written for each element.
 * \param [in]  write_mpirank If true, the mpirank is written for each element .
 * \param [in]  write_level If true, the refinement level is written for each element.
 * \param [in]  write_element_id If true, the global element id is written for each element.
 * \param [in]  write_ghosts If true, each process additionally writes its ghost elements.
 *                           For ghost element the treeid is -1.
 * \param [in]  num_data  Number of user defined double valued data fields to write.
 * \param [in]  data      Array of t8_vtk_data_field_t of length \a num_data
 *                        providing the used defined per element data.
 *                        If scalar and vector fields are used, all scalar fields
 *                        must come first in the array.
 * \return  True if succesful, false if not (process local).
 */
int                 t8_forest_vtk_write_file (t8_forest_t forest,
                                              const char *fileprefix,
                                              int write_treeid,
                                              int write_mpirank,
                                              int write_level,
                                              int write_element_id,
                                              int write_ghosts,
                                              int num_data,
                                              t8_vtk_data_field_t * data);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_VTK_H */
