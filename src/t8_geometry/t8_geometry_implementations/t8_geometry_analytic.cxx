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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.hxx>

t8_geometry_analytic::t8_geometry_analytic (int dim, std::string name, t8_geom_analytic_fn analytical,
                                            t8_geom_analytic_jacobian_fn jacobian_in,
                                            t8_geom_load_tree_data_fn load_tree_data_in,
                                            t8_geom_tree_negative_volume_fn tree_negative_volume_in,
                                            const void *user_data_in)
  : t8_geometry (dim, name + "_" + std::to_string (dim))
{
  analytical_function = analytical;
  jacobian = jacobian_in;
  load_tree_data = load_tree_data_in;
  tree_negative_volume = tree_negative_volume_in;
  user_data = user_data_in;
}

t8_geometry_analytic::t8_geometry_analytic (int dim, std::string name)
  : t8_geometry (dim, name + "_" + std::to_string (dim))
{
  analytical_function = NULL;
  jacobian = NULL;
  load_tree_data = NULL;
  tree_negative_volume = NULL;
  user_data = NULL;
}

void
t8_geometry_analytic::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                        const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (analytical_function != NULL);
  analytical_function (cmesh, gtreeid, ref_coords, num_coords, out_coords, tree_data, user_data);
}

void
t8_geometry_analytic::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                 const size_t num_coords, double *jacobian_out) const
{
  T8_ASSERT (jacobian != NULL);
  jacobian (cmesh, gtreeid, ref_coords, num_coords, jacobian_out, tree_data, user_data);
}

void
t8_geometry_analytic::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  if (load_tree_data != NULL) {
    /* Load tree data if a loading function was provided. */
    load_tree_data (cmesh, gtreeid, &tree_data);
  }
  else {
    /* Otherwise it is NULL. */
    tree_data = NULL;
  }
}

bool
t8_geometry_analytic::t8_geom_tree_negative_volume () const
{
  if (tree_negative_volume != NULL) {
    /* Tree negative volume if a loading function was provided. */
    return tree_negative_volume ();
  }
  else {
    return false;
  }
}

T8_EXTERN_C_BEGIN ();

void
t8_geometry_analytic_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);

  delete *geom;
  *geom = NULL;
}

t8_geometry_c *
t8_geometry_analytic_new (int dim, const char *name, t8_geom_analytic_fn analytical,
                          t8_geom_analytic_jacobian_fn jacobian, t8_geom_load_tree_data_fn load_tree_data,
                          t8_geom_tree_negative_volume_fn tree_negative_volume, const void *user_data)
{
  t8_geometry_analytic *geom
    = new t8_geometry_analytic (dim, name, analytical, jacobian, load_tree_data, tree_negative_volume, user_data);
  return (t8_geometry_c *) geom;
}

void
t8_geom_load_tree_data_vertices (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const void **user_data)
{
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  *user_data = t8_cmesh_get_tree_vertices (cmesh, ltreeid);
}

T8_EXTERN_C_END ();
