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

#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#include <t8_geometry/t8_geometry_helpers.h>

/**
 * Check that the two points of the geometry are ordered correctly, that is
 * p1_x <= p2_x, p1_y <= p2_y and p1_z <= p2_z
 * 
 * \param[in] tree_vertices The vertices of a tree
 * \return true if the points are ordered correctly
 * \return false ow
 */
#if T8_ENABLE_DEBUG
static bool
correct_point_order (const double tree_vertices[6])
{
  return tree_vertices[0] <= tree_vertices[3] && tree_vertices[1] <= tree_vertices[4]
         && tree_vertices[2] <= tree_vertices[5];
}
#endif

t8_geometry_linear_axis_aligned::t8_geometry_linear_axis_aligned ()
  : t8_geometry_with_vertices ("t8_geom_linear_axis_aligned")
{
}

t8_geometry_linear_axis_aligned::~t8_geometry_linear_axis_aligned ()
{
}

void
t8_geometry_linear_axis_aligned::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                   const size_t num_coords, double *out_coords) const
{
  T8_ASSERT (correct_point_order (active_tree_vertices));
  t8_geom_compute_linear_axis_aligned_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords,
                                                out_coords);
}

void
t8_geometry_linear_axis_aligned::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                                                            const double *ref_coords, const size_t num_coords,
                                                            double *jacobian) const
{
  SC_ABORT ("Not implemented.");
}

void
t8_geometry_linear_axis_aligned::t8_geom_point_batch_inside_element (t8_forest_t forest, t8_locidx_t ltreeid,
                                                                     const t8_element_t *element, const double *points,
                                                                     const int num_points, int *is_inside,
                                                                     const double tolerance) const
{
  double v_min[3];
  double v_max[3];

  const t8_eclass_t elem_class = t8_forest_get_tree_class (forest, ltreeid);
  T8_ASSERT (elem_class == T8_ECLASS_LINE || elem_class == T8_ECLASS_HEX || elem_class == T8_ECLASS_QUAD);

  const int max_corner_id = (elem_class == T8_ECLASS_LINE) ? 1 : (elem_class == T8_ECLASS_QUAD) ? 3 : 7;

  /*Geometry is fully described by v_min and v_max*/
  t8_forest_element_coordinate (forest, ltreeid, element, 0, v_min);
  t8_forest_element_coordinate (forest, ltreeid, element, max_corner_id, v_max);
#if T8_ENABLE_DEBUG
  const double coords[6] = { v_min[0], v_min[1], v_min[2], v_max[0], v_max[1], v_max[2] };
  T8_ASSERT (correct_point_order (coords));
#endif

  for (int ipoint = 0; ipoint < num_points; ipoint++) {
    /* A point is inside if it is inbetween the x/y/z-coordinates of v_min and v_max */
    /* check x-coordinate */
    /* check y-coordinate */
    /* check z-coordinate */
    is_inside[ipoint]
      = v_min[0] - tolerance <= points[ipoint * 3] && points[ipoint * 3] <= v_max[0] + tolerance
        && v_min[1] - tolerance <= points[ipoint * 3 + 1] && points[ipoint * 3 + 1] <= v_max[1] + tolerance
        && v_min[2] - tolerance <= points[ipoint * 3 + 2] && points[ipoint * 3 + 2] <= v_max[2] + tolerance;
  }
  return;
}

bool
t8_geometry_linear_axis_aligned::t8_geom_tree_negative_volume () const
{
  T8_ASSERT (correct_point_order (active_tree_vertices));
  return false;
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear_axis_aligned.h.
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_linear_axis_aligned_new ()
{
  t8_geometry_linear_axis_aligned *geom = new t8_geometry_linear_axis_aligned ();
  return (t8_geometry_c *) geom;
}

void
t8_geometry_linear_axis_aligned_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();
