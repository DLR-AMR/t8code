/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2023 the developers

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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_examples.hxx>
#include <t8_geometry/t8_geometry_helpers.h>
#include <t8_vec.h>

void
t8_geometry_squared_disk::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                            const size_t num_coords, double *out_coords) const
{
  if (num_coords != 1)
    SC_ABORT ("Error: Batch computation of geometry not yet supported.");

  double n[3]; /* Normal vector. */
  double r[3]; /* Radial vector. */
  double s[3]; /* Radial vector for the corrected coordinates. */
  double p[3]; /* Vector on the plane resp. quad. */

  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  double *tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  /* Center square. */
  if (gtreeid == 0) {

    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      size_t offset = 3 * i_coord;

      t8_geom_linear_interpolation (ref_coords + offset, tree_vertices, 3, 2, p);

      out_coords[offset + 0] = p[0];
      out_coords[offset + 1] = p[1];
      out_coords[offset + 2] = 0.0;
    }

    return;
  }

  /* Four squares framing the central one. */
  {
    const double center_ref[3] = { 0.5, 0.5, 0.0 };
    t8_geom_linear_interpolation (center_ref, tree_vertices, 3, 2, n);

    /* Normalize vector `n`. */
    const double norm = sqrt (n[0] * n[0] + n[1] * n[1]);
    n[0] = n[0] / norm;
    n[1] = n[1] / norm;
  }

  {
    /* Radial vector parallel to one of the tilted edges of the quad. */
    r[0] = tree_vertices[0];
    r[1] = tree_vertices[1];

    /* Normalize vector `r`. */
    const double norm = sqrt (r[0] * r[0] + r[1] * r[1]);
    r[0] = r[0] / norm;
    r[1] = r[1] / norm;
  }

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    size_t offset = 3 * i_coord;

    const double x = ref_coords[offset + 0];
    const double y = ref_coords[offset + 1];

    {
      double corr_ref_coords[3];

      /* Correction in order to rectify elements near the corners. */
      corr_ref_coords[0] = tan (0.5 * M_PI * (x - 0.5)) * 0.5 + 0.5;
      corr_ref_coords[1] = y;
      corr_ref_coords[2] = 0.0;

      /* Compute and normalize vector `s`. */
      t8_geom_linear_interpolation (corr_ref_coords, tree_vertices, 3, 2, s);

      const double norm = sqrt (s[0] * s[0] + s[1] * s[1]);
      s[0] = s[0] / norm;
      s[1] = s[1] / norm;
    }

    t8_geom_linear_interpolation (ref_coords + offset, tree_vertices, 3, 2, p);

    /* Compute intersection of line with a plane. */
    const double out_radius = (p[0] * n[0] + p[1] * n[1]) / (r[0] * n[0] + r[1] * n[1]);

    const double blend = y * out_radius; /* y \in [0,1] */
    const double dnelb = 1.0 - y;

    out_coords[offset + 0] = dnelb * p[0] + blend * s[0];
    out_coords[offset + 1] = dnelb * p[1] + blend * s[1];
    out_coords[offset + 2] = 0.0;
  }
}

T8_EXTERN_C_BEGIN ();

void
t8_geometry_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);

  delete *geom;
  *geom = NULL;
}

/* Satisfy the C interface from t8_geometry_linear.h. */
t8_geometry_c *
t8_geometry_squared_disk_new ()
{
  t8_geometry_squared_disk *geom = new t8_geometry_squared_disk ();
  return (t8_geometry_c *) geom;
}

T8_EXTERN_C_END ();
