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
#ifndef T8_VIS_H
#define T8_VIS_H

#include <t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8.h>

#if T8_WITH_VTK
#if T8_ENABLE_MPI
#include <vtkMultiProcessController.h>
#endif
#include <vtkUnstructuredGrid.h>

T8_EXTERN_C_BEGIN ();

typedef struct t8_interactive_vis
{
  /* Flag, if we have already read the Data from a file
   * non-zero, if the data has been read */
  int                 data_has_been_read;

  /* The refinement level of the forest. */
  int                 refinement_lvl;

  /* cmesh representing the data */
  t8_forest_t         forest;

  /* Name of the file to read. */
  char               *filepath;

  vtkUnstructuredGrid *vtkGrid;
  /* The communicator used. */
  sc_MPI_Comm         comm;
} t8_interactive_vis_t;

/**
 * Initalize an interactive visualization handler.
 * 
 * \param [in, out] pvis_handler The visualization handler to initialize.
 */
void                t8_interactive_vis_init (t8_interactive_vis_t **
                                             pvis_handler);

/**
 * Set the the filepath to the file to read the data from.
 * 
 * \param[in, out] vis_handler An initialized vis_handler.
 * \param[in] filepath    Path to the file.
 * \return int 
 */
void                t8_interactive_vis_set_filepath (t8_interactive_vis_t *
                                                     vis_handler,
                                                     const char *filepath);

/**
 * Set the pointer to the unstructured Grid 
 * 
 * \param[in, out] vis_handler An initialized vis_handler.
 * \param[in] grid             A vtkSmartPointer to an unstructuredGrid. 
 */
void
 
 
 
 
 
 
 
 t8_interactive_vis_set_vtkGrid (t8_interactive_vis_t * vis_handler,
                                 vtkSmartPointer < vtkUnstructuredGrid >
                                 grid);

/**
 * Set a uniform refinment of the forest.  
 * 
 * \param[in, out] vis_hanlder       An initialized vis_handler.
 * \param[in] level             The level we want to refine the forest to.
 */
void                t8_interactive_vis_set_refinement (t8_interactive_vis_t *
                                                       vis_hanlder,
                                                       const int level);

/**
 * Set the MPI communicator to use
 * 
 * \param[in, out] vis_handler  An initialized vis_handler. Its communicator will be set to \a comm
 * \param[in]      comm         An MPI-Communicator. 
 */
void                t8_interactive_vis_set_MPI_comm (t8_interactive_vis_t *
                                                     vis_handler,
                                                     sc_MPI_Comm comm);

/**
 * Update the vtkGrid of the vis_handler. If the Data has not been read already, we read
 * the Data from a vtk-file. 
 * If the Data has been read, we will do something else in the future (update the resolution for example)
 * 
 * \param vis_handler Initialized handler to update
 */
void                t8_interactive_vis_update_vtkGrid (t8_interactive_vis *
                                                       vis_handler);
/**
 * Destroy an interactive visualization handler. 
 * \param [in, out] pvis_hanlder The handler to destroy.
 * 
 */
void                t8_interactive_vis_destroy (t8_interactive_vis_t **
                                                pvis_handler);

T8_EXTERN_C_END ();

#endif /* T8_WITH_VTK */

#endif /* !T8_VIS_H */
