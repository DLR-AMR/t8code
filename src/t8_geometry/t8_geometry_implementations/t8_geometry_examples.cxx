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

/**
 * Map the faces of an oktaeder to a spherical surface.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_triangulated_spherical_surface::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                                                              const double *ref_coords, const size_t num_coords,
                                                              double *out_coords) const
{
  t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  double *tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  /* We average over the three corners of the triangle. */
  const double avg_factor = 1.0 / 3.0;

  /* Radius of the sphere scaled by the average factor. */
  const double radius = t8_vec_norm (tree_vertices) * avg_factor;

  /* The next three code blocks straighten out the elements near the triangle
   * corners by averaging the rectification with all three corners. */

  /* First triangle corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = tree_vertices[0];
    u[1] = tree_vertices[1];
    u[2] = tree_vertices[2];

    v[0] = tree_vertices[3 + 0] - u[0];
    v[1] = tree_vertices[3 + 1] - u[1];
    v[2] = tree_vertices[3 + 2] - u[2];

    w[0] = tree_vertices[6 + 0] - u[0];
    w[1] = tree_vertices[6 + 1] - u[1];
    w[2] = tree_vertices[6 + 2] - u[2];

    /* Reference coordinates from this particular triangle corner. */
    const double u_ref[3] = { 0.0, 0.0, 0.0 };
    const double v_ref[3] = { 1.0, 0.0, 0.0 };
    const double w_ref[3] = { -1.0, 1.0, 0.0 };

    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;

      const double x = ref_coords[offset + 0];
      const double y = ref_coords[offset + 1];

      /* Compute local triangle coordinate. */
      const double vv = u_ref[0] + x * v_ref[0] + y * w_ref[0];
      const double ww = u_ref[1] + x * v_ref[1] + y * w_ref[1];

      /* Correction in order to rectify elements near the corners. */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double ray[3]; /* Ray vector pinning through the triangle at reference coordinates. */
      ray[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      ray[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      ray[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      t8_vec_normalize (ray);

      out_coords[offset + 0] = radius * ray[0];
      out_coords[offset + 1] = radius * ray[1];
      out_coords[offset + 2] = radius * ray[2];
    }
  }

  /* Second triangle corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = tree_vertices[6 + 0];
    u[1] = tree_vertices[6 + 1];
    u[2] = tree_vertices[6 + 2];

    v[0] = tree_vertices[0 + 0] - u[0];
    v[1] = tree_vertices[0 + 1] - u[1];
    v[2] = tree_vertices[0 + 2] - u[2];

    w[0] = tree_vertices[3 + 0] - u[0];
    w[1] = tree_vertices[3 + 1] - u[1];
    w[2] = tree_vertices[3 + 2] - u[2];

    /* Reference coordinates from this particular triangle corner. */
    const double u_ref[3] = { 1.0, 0.0, 0.0 };
    const double v_ref[3] = { -1.0, 1.0, 0.0 };
    const double w_ref[3] = { 0.0, -1.0, 0.0 };

    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;

      const double x = ref_coords[offset + 0];
      const double y = ref_coords[offset + 1];

      /* Compute local triangle coordinate. */
      const double vv = u_ref[0] + x * v_ref[0] + y * w_ref[0];
      const double ww = u_ref[1] + x * v_ref[1] + y * w_ref[1];

      /* Correction in order to rectify elements near the corners. */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double ray[3]; /* Ray vector pinning through the triangle at reference coordinates. */
      ray[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      ray[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      ray[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      t8_vec_normalize (ray);

      out_coords[offset + 0] = out_coords[offset + 0] + radius * ray[0];
      out_coords[offset + 1] = out_coords[offset + 1] + radius * ray[1];
      out_coords[offset + 2] = out_coords[offset + 2] + radius * ray[2];
    }
  }

  /* Third triangle corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = tree_vertices[3 + 0];
    u[1] = tree_vertices[3 + 1];
    u[2] = tree_vertices[3 + 2];

    v[0] = tree_vertices[6 + 0] - u[0];
    v[1] = tree_vertices[6 + 1] - u[1];
    v[2] = tree_vertices[6 + 2] - u[2];

    w[0] = tree_vertices[0 + 0] - u[0];
    w[1] = tree_vertices[0 + 1] - u[1];
    w[2] = tree_vertices[0 + 2] - u[2];

    /* Reference coordinates from this particular triangle corner. */
    const double u_ref[3] = { 0.0, 1.0, 0.0 };
    const double v_ref[3] = { 0.0, -1.0, 0.0 };
    const double w_ref[3] = { 1.0, 0.0, 0.0 };

    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;

      const double x = ref_coords[offset + 0];
      const double y = ref_coords[offset + 1];

      /* Compute local triangle coordinate. */
      const double vv = u_ref[0] + x * v_ref[0] + y * w_ref[0];
      const double ww = u_ref[1] + x * v_ref[1] + y * w_ref[1];

      /* Correction in order to rectify elements near the corners. */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double ray[3]; /* Ray vector pinning through the triangle at reference coordinates. */
      ray[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      ray[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      ray[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      t8_vec_normalize (ray);

      out_coords[offset + 0] = out_coords[offset + 0] + radius * ray[0];
      out_coords[offset + 1] = out_coords[offset + 1] + radius * ray[1];
      out_coords[offset + 2] = out_coords[offset + 2] + radius * ray[2];
    }
  }
}

/**
 * Map the faces of a unit cube to a spherical surface.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_quadrangulated_spherical_surface::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid,
                                                                const double *ref_coords, const size_t num_coords,
                                                                double *out_coords) const
{
  double n[3]; /* Normal vector. */
  double r[3]; /* Radial vector. */
  double p[3]; /* Vector on the plane. */

  const t8_locidx_t ltreeid = t8_cmesh_get_local_id (cmesh, gtreeid);
  const double *tree_vertices = t8_cmesh_get_tree_vertices (cmesh, ltreeid);

  t8_geom_linear_interpolation (t8_element_centroid_ref_coords[T8_ECLASS_QUAD], tree_vertices, 3, 2, n);
  t8_vec_normalize (n);

  r[0] = tree_vertices[0];
  r[1] = tree_vertices[1];
  r[2] = tree_vertices[2];

  t8_vec_normalize (r);

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset = 3 * i_coord;

    {
      double corr_ref_coords[3]; /* Corrected reference coordinates. */

      const double x = ref_coords[offset + 0];
      const double y = ref_coords[offset + 1];
      const double z = ref_coords[offset + 2];

      /* tldr: Correction in order to rectify elements near the corners. 
       * This is necessary, since due to the transformation from the unit cube
       * to the sphere elements near the face centers expand while near the
       * corners they shrink. Following correction alleviates this.
       */
      corr_ref_coords[0] = tan (0.5 * M_PI * (x - 0.5)) * 0.5 + 0.5;
      corr_ref_coords[1] = tan (0.5 * M_PI * (y - 0.5)) * 0.5 + 0.5;
      corr_ref_coords[2] = z;

      t8_geom_linear_interpolation (corr_ref_coords, tree_vertices, 3, 2, p);
    }

    const double R = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / (r[0] * n[0] + r[1] * n[1] + r[2] * n[2]);

    t8_vec_normalize (p);

    out_coords[offset + 0] = R * p[0];
    out_coords[offset + 1] = R * p[1];
    out_coords[offset + 2] = R * p[2];
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

t8_geometry_c *
t8_geometry_triangulated_spherical_surface_new ()
{
  t8_geometry_triangulated_spherical_surface *geom = new t8_geometry_triangulated_spherical_surface ();
  return (t8_geometry_c *) geom;
}

t8_geometry_c *
t8_geometry_quadrangulated_spherical_surface_new ()
{
  t8_geometry_quadrangulated_spherical_surface *geom = new t8_geometry_quadrangulated_spherical_surface ();
  return (t8_geometry_c *) geom;
}

T8_EXTERN_C_END ();
