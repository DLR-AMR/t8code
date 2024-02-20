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
#include <algorithm>

t8_geometry_handler::t8_geometry_handler (): active_geometry (nullptr), active_tree (-1)
{
}

t8_geometry_handler::~t8_geometry_handler ()
{
  /* Nothing to do */
}

template <typename geometry, typename... args>
t8_geometry &
t8_geometry_handler::register_geometry (args &&...args)
{
  std::unique_ptr<t8_geometry> geom = std::make_unique<geometry> (std::forward<Args> (args)...);
  const size_t hash = geom->t8_geom_get_hash ();
  if (registered_geometries.find (hash) == registered_geometries.end ()) {
    registered_geometries.emplace (hash, geom);
  }
  if (registered_geometries.size () == 1) {
    active_geometry = registered_geometries.at (hash).get ();
  }
  return *registered_geometries.at (hash).get ();
}

t8_geometry &
t8_geometry_handler::register_geometry (t8_geometry &geom)
{
  const size_t hash = geom.t8_geom_get_hash ();
  if (registered_geometries.find (hash) == registered_geometries.end ()) {
    registered_geometries.emplace (hash, std::make_unique<t8_geometry> (std::move (geom)));
  }
  if (registered_geometries.size () == 1) {
    active_geometry = registered_geometries.at (hash).get ();
  }
  return *registered_geometries.at (hash).get ();
}

inline t8_geometry *
t8_geometry_handler::get_geometry (const std::string &name)
{
  const size_t hash = std::hash<std::string> {}(name);
  return t8_geometry_handler::get_geometry (hash);
}

inline t8_geometry *
t8_geometry_handler::get_geometry (const size_t hash)
{
  auto found = registered_geometries.find (hash);
  if (found != registered_geometries.end ()) {
    return found->second.get ();
  }
  return nullptr;
}

inline t8_geometry *
t8_geometry_handler::get_unique_geometry ()
{
  T8_ASSERT (registered_geometries.size () == 1);
  return active_geometry;
}

inline t8_geometry *
t8_geometry_handler::get_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  update_tree (cmesh, gtreeid);
  return active_geometry;
}

inline void
t8_geometry_handler::evaluate_tree_geometry (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                             const size_t num_coords, double *out_coords)
{
  update_tree (cmesh, gtreeid);
  active_geometry->t8_geom_evaluate (cmesh, gtreeid, ref_coords, num_coords, out_coords);
}

inline void
t8_geometry_handler::evaluate_tree_geometry_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                      const size_t num_coords, double *out_coords)
{
  update_tree (cmesh, gtreeid);
  active_geometry->t8_geom_evaluate_jacobian (cmesh, gtreeid, ref_coords, num_coords, out_coords);
}

inline t8_geometry_type_t
t8_geometry_handler::get_tree_geometry_type (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  update_tree (cmesh, gtreeid);
  return active_geometry->t8_geom_get_type ();
}

void
t8_geometry_handler::update_tree (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  T8_ASSERT (0 <= gtreeid && gtreeid < t8_cmesh_get_num_trees (cmesh));
  if (active_tree != gtreeid) {
    const int num_geoms = get_num_geometries ();
    /* This tree is not the active tree. We need to update the 
     * active tree, its geometry and its data. */
    /* Set the new tree as active. */
    active_tree = gtreeid;
    if (num_geoms > 1) {
      /* Find and load the geometry of that tree. 
       * Only necessary if we have more than one geometry. */
      const std::string geom_name (t8_cmesh_get_tree_geom_name (cmesh, gtreeid));
      active_geometry = get_geometry (geom_name);
      SC_CHECK_ABORTF (active_geometry != nullptr,
                       "Could not find geometry with name %s or tree %d has no registered geometry.",
                       geom_name.c_str (), gtreeid);
    }
    /* Get the user data for this geometry and this tree. */
    active_geometry->t8_geom_load_tree_data (cmesh, gtreeid);
  }
  T8_ASSERT (active_geometry != nullptr);
}