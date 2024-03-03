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

#include <t8_geometry/t8_geometry_handler.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_geometry.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_cad.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>

#include <algorithm>
#include <memory>

void
t8_geometry_handler::register_geometry (t8_geometry_c **geom)
{
  std::unique_ptr<t8_geometry> geom_ptr = std::unique_ptr<t8_geometry> (std::move (*geom));
  *geom = add_geometry<t8_geometry> (std::move (geom_ptr));
}

void
t8_geometry_handler::update_tree (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
#if T8_ENABLE_DEBUG
static int num_recursive_invocations = 0;
  /*******************************************************/
  /* We ran into a bug, where t8_cmesh_is_committed called update_tree for a different tree,
   * so that at the end of this function the tree with gtreeid actually was not loaded.
   * The do_not_update_tree variable lives beyond of this function call and is used to
   * alarm if the update_tree function is being called recursively. */
  T8_ASSERT(t8_cmesh_is_committed(cmesh));
  static int do_not_update_tree = false;
  static bool is_checking = false; /* This variable lives beyond of the scope of this function and prevents t8_cmesh_is_committed from executing */
  t8_debugf ("update_tree called %i times recursively\n", num_recursive_invocations++);
  t8_debugf ("is_checking: %i\n", is_checking);
  t8_debugf ("do_not_update_tree: %i\n", do_not_update_tree);
  T8_ASSERT (!do_not_update_tree); /* If this assertion is triggers, the update_tree function is being called recursively. */
  do_not_update_tree = true;
  const bool former_state_of_is_checking = is_checking;
  is_checking = true;
  /*******************************************************/
#endif /* T8_ENABLE_DEBUG */
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_cmesh_get_num_trees (cmesh));
  const int num_geoms = get_num_geometries ();
  SC_CHECK_ABORTF (num_geoms > 0,
                   "The geometry of the tree could not be loaded, because no geometries were registered.");
  T8_ASSERT (active_geometry != nullptr);
  if (active_tree != gtreeid) {
    /* This tree is not the active tree. We need to update the 
     * active tree, its geometry and its data. */
    /* Set the new tree as active. */
    active_tree = gtreeid;
    if (num_geoms > 1) {
      /* Find and load the geometry of that tree. 
       * Only necessary if we have more than one geometry. */
      const size_t geom_hash = t8_cmesh_get_tree_geom_hash (cmesh, gtreeid);
      active_geometry = get_geometry (geom_hash);
      SC_CHECK_ABORTF (active_geometry != nullptr,
                       "Could not find geometry with hash %zu or tree %ld has no registered geometry.", geom_hash,
                       gtreeid);
    }
    /* Get the user data for this geometry and this tree. */
    active_geometry->t8_geom_load_tree_data (cmesh, gtreeid);
  }
#if T8_ENABLE_DEBUG
  is_checking = former_state_of_is_checking;
  do_not_update_tree = false;
  num_recursive_invocations = 0;
#endif /* T8_ENABLE_DEBUG */
}