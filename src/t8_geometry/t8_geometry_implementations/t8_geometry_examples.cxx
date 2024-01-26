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

  /* Center square. */
  if (gtreeid == 0) {

    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      size_t offset = 3 * i_coord;

      t8_geom_linear_interpolation (ref_coords + offset, active_tree_vertices, 3, 2, p);

      out_coords[offset + 0] = p[0];
      out_coords[offset + 1] = p[1];
      out_coords[offset + 2] = 0.0;
    }

    return;
  }

  /* Four squares framing the central one. */
  {
    const double center_ref[3] = { 0.5, 0.5, 0.0 };
    t8_geom_linear_interpolation (center_ref, active_tree_vertices, 3, 2, n);

    /* Normalize vector `n`. */
    const double norm = sqrt (n[0] * n[0] + n[1] * n[1]);
    n[0] = n[0] / norm;
    n[1] = n[1] / norm;
  }

  {
    /* Radial vector parallel to one of the tilted edges of the quad. */
    r[0] = active_tree_vertices[0];
    r[1] = active_tree_vertices[1];

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
      t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, 2, s);

      const double norm = sqrt (s[0] * s[0] + s[1] * s[1]);
      s[0] = s[0] / norm;
      s[1] = s[1] / norm;
    }

    t8_geom_linear_interpolation (ref_coords + offset, active_tree_vertices, 3, 2, p);

    /* Compute intersection of line with a plane. */
    const double out_radius = (p[0] * n[0] + p[1] * n[1]) / (r[0] * n[0] + r[1] * n[1]);

    const double blend = y * out_radius; /* y \in [0,1] */
    const double dnelb = 1.0 - y;

    out_coords[offset + 0] = dnelb * p[0] + blend * s[0];
    out_coords[offset + 1] = dnelb * p[1] + blend * s[1];
    out_coords[offset + 2] = 0.0;
  }
}

static inline void
t8_geom_evaluate_sphere_tri_prism (const double *active_tree_vertices, const t8_eclass_t eclass,
                                   const double *ref_coords, const size_t num_coords, double *out_coords)
{
  double n[3]; /* Normal vector of the current triangle. For prisms along z-axis in reference space. */
  t8_vec_tri_normal (active_tree_vertices, active_tree_vertices + 3, active_tree_vertices + 6, n);
  t8_vec_normalize (n);

  double r[3]; /* Radial vector through one of triangle/prism's corners. */
  r[0] = active_tree_vertices[0];
  r[1] = active_tree_vertices[1];
  r[2] = active_tree_vertices[2];
  t8_vec_normalize (r);

  /* The next three code blocks straighten out the elements near the triangle
   * corners by averaging the rectification with all three corners. */

  /* With this factor we compute the intersection of `r` and `p` (see further
   * below) and average over the three corners of the triangle. */
  const double factor = 1.0 / (r[0] * n[0] + r[1] * n[1] + r[2] * n[2]) / 3.0;

  /* First triangle/prism corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = active_tree_vertices[0];
    u[1] = active_tree_vertices[1];
    u[2] = active_tree_vertices[2];

    v[0] = active_tree_vertices[3 + 0] - u[0];
    v[1] = active_tree_vertices[3 + 1] - u[1];
    v[2] = active_tree_vertices[3 + 2] - u[2];

    w[0] = active_tree_vertices[6 + 0] - u[0];
    w[1] = active_tree_vertices[6 + 1] - u[1];
    w[2] = active_tree_vertices[6 + 2] - u[2];

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

      /* tldr: Correction in order to rectify elements near the corners. This
       * is necessary, since due to the transformation from the cmesh triangle to the
       * sphere elements near the face centers expand while near the corners they
       * shrink. Following correction alleviates this.
       * TODO: This correction is not general and not optimal in all cases.
       *       Find a better one. (This is not a trivial task, though.)
       */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double p[3];
      p[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      p[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      p[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      const double radius = factor * (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / t8_vec_norm (p);

      out_coords[offset + 0] = radius * p[0];
      out_coords[offset + 1] = radius * p[1];
      out_coords[offset + 2] = radius * p[2];
    }
  }

  /* Second triangle/prism corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = active_tree_vertices[6 + 0];
    u[1] = active_tree_vertices[6 + 1];
    u[2] = active_tree_vertices[6 + 2];

    v[0] = active_tree_vertices[0 + 0] - u[0];
    v[1] = active_tree_vertices[0 + 1] - u[1];
    v[2] = active_tree_vertices[0 + 2] - u[2];

    w[0] = active_tree_vertices[3 + 0] - u[0];
    w[1] = active_tree_vertices[3 + 1] - u[1];
    w[2] = active_tree_vertices[3 + 2] - u[2];

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

      /* tldr: Correction in order to rectify elements near the corners. This
       * is necessary, since due to the transformation from the cmesh triangle to the
       * sphere elements near the face centers expand while near the corners they
       * shrink. Following correction alleviates this.
       * TODO: This correction is not general and not optimal in all cases.
       *       Find a better one. (This is not a trivial task, though.)
       */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double p[3];
      p[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      p[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      p[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      const double norm = t8_vec_norm (p);
      const double R = factor * (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / norm;

      /* Note, in `R` there already is the avg. factor `1/3` included. */
      out_coords[offset + 0] = out_coords[offset + 0] + R * p[0];
      out_coords[offset + 1] = out_coords[offset + 1] + R * p[1];
      out_coords[offset + 2] = out_coords[offset + 2] + R * p[2];
    }
  }

  /* Third triangle/prism corner. */
  {
    double u[3]; /* Position vector. */
    double v[3]; /* First triangle side. */
    double w[3]; /* Second triangle side. */

    u[0] = active_tree_vertices[3 + 0];
    u[1] = active_tree_vertices[3 + 1];
    u[2] = active_tree_vertices[3 + 2];

    v[0] = active_tree_vertices[6 + 0] - u[0];
    v[1] = active_tree_vertices[6 + 1] - u[1];
    v[2] = active_tree_vertices[6 + 2] - u[2];

    w[0] = active_tree_vertices[0 + 0] - u[0];
    w[1] = active_tree_vertices[0 + 1] - u[1];
    w[2] = active_tree_vertices[0 + 2] - u[2];

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

      /* tldr: Correction in order to rectify elements near the corners. This
       * is necessary, since due to the transformation from the cmesh triangle to the
       * sphere elements near the face centers expand while near the corners they
       * shrink. Following correction alleviates this.
       * TODO: This correction is not general and not optimal in all cases.
       *       Find a better one. (This is not a trivial task, though.)
       */
      const double vv_corr = tan (0.5 * M_PI * (vv - 0.5)) * 0.5 + 0.5;
      const double ww_corr = tan (0.5 * M_PI * (ww - 0.5)) * 0.5 + 0.5;

      /* Compute and apply the corrected mapping. */
      double p[3];
      p[0] = u[0] + vv_corr * v[0] + ww_corr * w[0];
      p[1] = u[1] + vv_corr * v[1] + ww_corr * w[1];
      p[2] = u[2] + vv_corr * v[2] + ww_corr * w[2];

      const double norm = t8_vec_norm (p);
      const double R = factor * (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / norm;

      /* Note, in `R` there already is the avg. factor `1/3` included. */
      out_coords[offset + 0] = out_coords[offset + 0] + R * p[0];
      out_coords[offset + 1] = out_coords[offset + 1] + R * p[1];
      out_coords[offset + 2] = out_coords[offset + 2] + R * p[2];
    }
  }

  /* For triangles we are done. */
  if (eclass == T8_ECLASS_TRIANGLE)
    return;

  /* With this factor we compute the intersection of `r` and `p` (see further below). */
  const double denominator = 1.0 / (r[0] * n[0] + r[1] * n[1] + r[2] * n[2]);

  /* We loop again over all points and correct the length in radial direction to pad the shell thickness. */
  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset = 3 * i_coord;

    double p[3];
    t8_geom_compute_linear_geometry (T8_ECLASS_PRISM, active_tree_vertices, ref_coords + offset, 1, p);
    const double current_radius = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) * denominator;

    t8_vec_normalize (out_coords + offset);
    out_coords[offset + 0] = out_coords[offset + 0] * current_radius;
    out_coords[offset + 1] = out_coords[offset + 1] * current_radius;
    out_coords[offset + 2] = out_coords[offset + 2] * current_radius;
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
  t8_geom_evaluate_sphere_tri_prism (active_tree_vertices, T8_ECLASS_TRIANGLE, ref_coords, num_coords, out_coords);
}

/**
 * Map the prismed faces of an oktaeder to a spherical shell.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_prismed_spherical_shell::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                       const size_t num_coords, double *out_coords) const

{
  t8_geom_evaluate_sphere_tri_prism (active_tree_vertices, T8_ECLASS_PRISM, ref_coords, num_coords, out_coords);
}

static inline void
t8_geom_evaluate_sphere_quad_hex (const double *active_tree_vertices, const int ndims, const double *ref_coords,
                                  const size_t num_coords, double *out_coords)
{
  double n[3]; /* Normal vector. */
  double r[3]; /* Radial vector. */
  double p[3]; /* Vector on the plane. */

  t8_geom_linear_interpolation (t8_element_centroid_ref_coords[T8_ECLASS_QUAD], active_tree_vertices, 3, 2, n);
  t8_vec_normalize (n);

  r[0] = active_tree_vertices[0];
  r[1] = active_tree_vertices[1];
  r[2] = active_tree_vertices[2];

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

      t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, ndims, p);
    }

    const double radius = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2]) / (r[0] * n[0] + r[1] * n[1] + r[2] * n[2]);

    t8_vec_normalize (p);

    out_coords[offset + 0] = radius * p[0];
    out_coords[offset + 1] = radius * p[1];
    out_coords[offset + 2] = radius * p[2];
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
  t8_geom_evaluate_sphere_quad_hex (active_tree_vertices, 2, ref_coords, num_coords, out_coords);
}

/**
 * Maps six hexaeders arranged into cube to a spherical shell.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_cubed_spherical_shell::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                     const size_t num_coords, double *out_coords) const
{
  t8_geom_evaluate_sphere_quad_hex (active_tree_vertices, 3, ref_coords, num_coords, out_coords);
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

t8_geometry_c *
t8_geometry_cubed_spherical_shell_new ()
{
  t8_geometry_cubed_spherical_shell *geom = new t8_geometry_cubed_spherical_shell ();
  return (t8_geometry_c *) geom;
}

t8_geometry_c *
t8_geometry_prismed_spherical_shell_new ()
{
  t8_geometry_prismed_spherical_shell *geom = new t8_geometry_prismed_spherical_shell ();
  return (t8_geometry_c *) geom;
}

T8_EXTERN_C_END ();
