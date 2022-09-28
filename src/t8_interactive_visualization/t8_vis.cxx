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

#include <t8_interactive_visualization/t8_vis.hxx>
#include <t8_cmesh/t8_cmesh_vtk_helper.hxx>

#if T8_WITH_VTK
#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>
#include <vtkMPI.h>
#endif

void
t8_interactive_vis_init (t8_interactive_vis_t ** pvis_handler)
{
  T8_ASSERT (pvis_handler != NULL);
  t8_interactive_vis_t *vis_handler = T8_ALLOC (t8_interactive_vis_t, 1);
  vis_handler->cmesh = NULL;
  vis_handler->data_has_been_read = 0;
  /* Is this a good value to intitialize the communicator? */
  vis_handler->comm = sc_MPI_COMM_NULL;
#if T8_WITH_VTK
  vis_handler->vtkGrid = NULL;
#endif
  *pvis_handler = vis_handler;
}

void
t8_interactive_vis_set_filenpath (t8_interactive_vis_t * vis_handler,
                                  char *filepath)
{
  strcpy (vis_handler->filepath, filepath);
}

void
t8_interactive_vis_set_MPI_comm (t8_interactive_vis_t * vis_handler,
                                 sc_MPI_Comm comm)
{
  vis_handler->comm = comm;
}

void
t8_interactive_vis_update_vtkGrid (t8_interactive_vis * vis_handler)
{
  if (!vis_handler->data_has_been_read) {
    /* TODO: Currently done twice, extrad t8_read_unstructured from t8_cmesh_read */
    t8_read_unstructured (vis_handler->filepath, vis_handler->vtkGrid);
    t8_cmesh_read_from_vtk_unstructured (vis_handler->filepath, 0, 0,
                                         vis_handler->comm);
  }
}

void
t8_interactive_vis_destroy (t8_interactive_vis_t ** pvis_handler)
{
  T8_ASSERT (pvis_handler != NULL);
  t8_interactive_vis_t *vis_handler = *pvis_handler;
  /* unreference/destroy the cmesh */
  t8_cmesh_unref (&vis_handler->cmesh);
  /* Free the memory */
  T8_FREE (vis_handler);
  /* Set the pointer to NULL */
  *pvis_handler = NULL;
}
