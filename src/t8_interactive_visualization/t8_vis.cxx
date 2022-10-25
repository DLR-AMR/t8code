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
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest.h>

#if T8_WITH_VTK
#include <t8_forest/t8_forest_vtk_helper.hxx>
#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>
#include <vtkMPI.h>

#if 0
int
t8_interactive_vis_adapt_callback (t8_forest_t forest,
                                   t8_forest_t forest_from,
                                   t8_locidx_t which_tree,
                                   t8_locidx_t lelement_id,
                                   t8_eclass_scheme_c *ts,
                                   const int is_family,
                                   const int num_elements,
                                   t8_element_t *elements[])
{
  if (ts->t8_element_level (elements[0]) < 1) {
    return 1;
  }
  else {
    return 0;
  }
}

void
t8_interactive_vis_update_vtkGrid (t8_interactive_vis_t * vis_handler)
{
  if (!vis_handler->data_has_been_read) {
    /* TODO: Currently done twice, extrad t8_read_unstructured from t8_cmesh_read */
    const int           successful_read =
      t8_read_unstructured (vis_handler->filepath, vis_handler->vtkGrid, 1, 0,
                            vis_handler->comm);
    t8_cmesh_t          cmesh;
    t8_cmesh_t          cmesh_in;
    t8_cmesh_init (&cmesh_in);
    t8_cmesh_init (&cmesh);
    if (successful_read) {
      t8_unstructured_to_cmesh (vis_handler->vtkGrid, 1, 0, cmesh_in,
                                vis_handler->comm);
      t8_cmesh_commit (cmesh_in, vis_handler->comm);
      t8_cmesh_set_derive (cmesh, cmesh_in);
      t8_cmesh_set_partition_uniform (cmesh, 0, t8_scheme_new_default_cxx ());
    }
    if (cmesh != NULL) {
      t8_cmesh_commit (cmesh, vis_handler->comm);
    }
    else {
      t8_global_errorf ("Could not commit cmesh.\n");
    }
    t8_forest_init (&vis_handler->forest);
    t8_forest_set_cmesh (vis_handler->forest, cmesh, vis_handler->comm);
    t8_forest_set_scheme (vis_handler->forest, t8_scheme_new_default_cxx ());
    t8_forest_commit (vis_handler->forest);
  }
  if (vis_handler->refinement_lvl > 0) {
    vis_handler->forest =
      t8_forest_new_adapt (vis_handler->forest,
                           t8_interactive_vis_adapt_callback, 0, 0, NULL);
  }
  /* TODO: Currently no data-writing. Enable writing of data. */
  t8_forest_to_vtkUnstructuredGrid (vis_handler->forest, vis_handler->vtkGrid,
                                    1, 1, 1, 1, 0, 0, NULL);
  t8_forest_write_vtk (vis_handler->forest, "vis_handler_");
}

void
t8_interactive_vis_destroy (t8_interactive_vis_t ** pvis_handler)
{
  T8_ASSERT (pvis_handler != NULL);
  t8_interactive_vis_t *vis_handler = *pvis_handler;
  /* unreference/destroy the cmesh */
  t8_forest_unref (&vis_handler->forest);
  T8_FREE (vis_handler->filepath);
  /* Free the memory */
  T8_FREE (vis_handler);
  /* Set the pointer to NULL */
  *pvis_handler = NULL;
  t8_debugf ("[D] destroyed the vis handler\n");
}
#endif
#endif
