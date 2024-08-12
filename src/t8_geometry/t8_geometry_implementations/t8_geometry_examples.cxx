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
t8_geometry_quadrangulated_disk::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                                   const size_t num_coords, double *out_coords) const
{
  double n[3] = { 0.0 }; /* Normal vector. */
  double r[3] = { 0.0 }; /* Radial vector. */
  double s[3] = { 0.0 }; /* Radial vector for the corrected coordinates. */
  double p[3] = { 0.0 }; /* Vector on the plane resp. quad. */

  /* Center quads. */
  if (gtreeid % 3 == 0) {
    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset_2d = 2 * i_coord;
      const size_t offset_3d = 3 * i_coord;
      t8_geom_linear_interpolation (ref_coords + offset_2d, active_tree_vertices, 3, 2, out_coords + offset_3d);
    }
    return;
  }

  /* Normal vector along one of the straight edges of the quad. */
  t8_vec_copy (active_tree_vertices, n);
  t8_vec_normalize (n);

  /* Radial vector parallel to one of the tilted edges of the quad. */
  t8_vec_copy (active_tree_vertices + 9, r);
  t8_vec_normalize (r);

  const double inv_denominator = 1.0 / t8_vec_dot (r, n);

  /* Radial reference coordinate index. */
  const int r_coord = ((gtreeid - 2) % 3 == 0) ? 0 : 1;
  /* Angular reference coordinate idex. */
  const int a_coord = ((gtreeid - 2) % 3 == 0) ? 1 : 0;

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset_2d = 2 * i_coord;
    const size_t offset_3d = 3 * i_coord;

    const double r_ref = ref_coords[offset_2d + r_coord];
    const double a_ref = ref_coords[offset_2d + a_coord];

    {
      double corr_ref_coords[3];

      /* Correction in order to rectify elements near the corners. */
      corr_ref_coords[r_coord] = r_ref;
      corr_ref_coords[a_coord] = tan (0.25 * M_PI * a_ref);
      corr_ref_coords[2] = 0.0;

      /* Compute and normalize vector `s`. */
      t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, 2, s);
      t8_vec_normalize (s);
    }

    /* Correction in order to rectify elements near the corners. */
    t8_geom_linear_interpolation (ref_coords + offset_2d, active_tree_vertices, 3, 2, p);

    /* Compute intersection of line with a plane. */
    const double out_radius = t8_vec_dot (p, n) * inv_denominator;

    /* Linear blend from flat to curved: `out_coords = (1.0 - r_ref)*p + r_ref_ * out_radius * s`. */
    t8_vec_axy (p, out_coords + offset_3d, 1.0 - r_ref);
    t8_vec_axpy (s, out_coords + offset_3d, r_ref * out_radius);
  }
}

static inline void
t8_geom_evaluate_sphere_tri_prism (const double *active_tree_vertices, const t8_eclass_t eclass,
                                   const double *ref_coords, const size_t num_coords, double *out_coords)
{
  // All elements are aligned such that the reference z-direction follows the
  // outward radial direction of the sphere. Hence the inner radius is equal to
  // the norm of the first positition vector of `active_tree_vertices`.
  const double inner_radius = t8_vec_norm (active_tree_vertices);

  t8_geom_compute_linear_geometry (eclass, active_tree_vertices, ref_coords, num_coords, out_coords);

  if (eclass == T8_ECLASS_TRIANGLE) {
    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;
      t8_vec_rescale (out_coords + offset, inner_radius);
    }
  }
  else {
    const size_t outer_vertex_offset = 3 * 3;
    const double shell_thickness = t8_vec_norm (active_tree_vertices + outer_vertex_offset) - inner_radius;
    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;
      const double z = ref_coords[offset + 2];
      t8_vec_rescale (out_coords + offset, inner_radius + z * shell_thickness);
    }
  }
}

/**
 * Map the faces of an octahedron to a spherical surface.
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
 * Map the prismed faces of an octahedron to a spherical shell.
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
  double position[3]; /* Position vector in the element. */

  /* All elements are aligned such that the face normal follows the
   * outward radial direction of the sphere. */
  const double radius = t8_vec_norm (active_tree_vertices);

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset_2d = 2 * i_coord;
    const size_t offset_3d = 3 * i_coord;

    double corr_ref_coords[3]; /* Corrected reference coordinates. */

    /* Shorthand for code readability. `ref_coords` go from 0 to 1. */
    const double x = ref_coords[offset_2d + 0];
    const double y = ref_coords[offset_2d + 1];

    /* tldr: Correction in order to rectify elements near the corners. 
     * This is necessary, since due to the transformation from the unit cube
     * to the sphere elements near the face centers expand while near the
     * corners they shrink. Following correction alleviates this.
     */
    corr_ref_coords[0] = tan (0.5 * M_PI * (x - 0.5)) * 0.5 + 0.5;
    corr_ref_coords[1] = tan (0.5 * M_PI * (y - 0.5)) * 0.5 + 0.5;
    corr_ref_coords[2] = 0;

    t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, 2, position);
    t8_vec_normalize (position);
    t8_vec_axy (position, out_coords + offset_3d, radius);
  }
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
  double position[3]; /* Position vector in the element. */

  /* All elements are aligned such that the reference z-direction follows the
   * outward radial direction of the sphere. Hence the element height is equal to
   * the shell thickness. */
  const double inner_radius = t8_vec_norm (active_tree_vertices);
  const double shell_thickness = t8_vec_norm (active_tree_vertices + 4 * 3) - inner_radius;

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset = 3 * i_coord;

    double corr_ref_coords[3]; /* Corrected reference coordinates. */

    /* Shorthand for code readability. `ref_coords` go from 0 to 1. */
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

    t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, 3, position);
    t8_vec_normalize (position);
    t8_vec_axy (position, out_coords + offset, inner_radius + z * shell_thickness);
  }
}

void
t8_geometry_cubed_sphere::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                            const size_t num_coords, double *out_coords) const
{
  double n[3] = { 0.0 }; /* Normal vector. */
  double r[3] = { 0.0 }; /* Radial vector. */
  double s[3] = { 0.0 }; /* Radial vector for the corrected coordinates. */
  double p[3] = { 0.0 }; /* Vector on the plane resp. quad. */

  /* Center hex. */
  if (gtreeid % 4 == 0) {
    for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
      const size_t offset = 3 * i_coord;
      t8_geom_linear_interpolation (ref_coords + offset, active_tree_vertices, 3, 3, out_coords + offset);
    }
    return;
  }

  t8_vec_copy (active_tree_vertices, n);
  t8_vec_normalize (n);

  t8_vec_copy (active_tree_vertices + 7 * 3, r);
  t8_vec_normalize (r);

  const double inv_denominator = 1.0 / t8_vec_dot (r, n);

  /* Radial reference coordinate index. */
  const int r_coord_lookup[4] = { -1, 1, 0, 2 };
  const int r_coord = r_coord_lookup[gtreeid % 4];
  /* Angular (theta) reference coordinate index. */
  const int t_coord_lookup[4] = { -1, 0, 1, 0 };
  const int t_coord = t_coord_lookup[gtreeid % 4];
  /* Angular (phi) reference coordinate index. */
  const int p_coord_lookup[4] = { -1, 2, 2, 1 };
  const int p_coord = p_coord_lookup[gtreeid % 4];

  for (size_t i_coord = 0; i_coord < num_coords; i_coord++) {
    const size_t offset = 3 * i_coord;

    const double r_ref = ref_coords[offset + r_coord]; /* radius */
    const double t_ref = ref_coords[offset + t_coord]; /* theta */
    const double p_ref = ref_coords[offset + p_coord]; /* phi */

    {
      double corr_ref_coords[3];

      /* Correction in order to rectify elements near the corners. Note, this
       * is probably not the most accurate correction but it does a decent enough job.
       */
      corr_ref_coords[r_coord] = r_ref;
      corr_ref_coords[t_coord] = tan (0.25 * M_PI * t_ref);
      corr_ref_coords[p_coord] = tan (0.25 * M_PI * p_ref);

      /* Compute and normalize vector `s`. */
      t8_geom_linear_interpolation (corr_ref_coords, active_tree_vertices, 3, 3, s);
      t8_vec_normalize (s);
    }

    t8_geom_linear_interpolation (ref_coords + offset, active_tree_vertices, 3, 3, p);

    /* Compute intersection of line with a plane. */
    const double out_radius = t8_vec_dot (p, n) * inv_denominator;

    /* Linear blend from flat to curved: `out_coords = (1.0 - r_ref)*p + r_ref_ * out_radius * s`. */
    t8_vec_axy (p, out_coords + offset, 1.0 - r_ref);
    t8_vec_axpy (s, out_coords + offset, r_ref * out_radius);
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
t8_geometry_quadrangulated_disk_new ()
{
  t8_geometry_quadrangulated_disk *geom = new t8_geometry_quadrangulated_disk ();
  return (t8_geometry_c *) geom;
}

t8_geometry_c *
t8_geometry_triangulated_spherical_surface_new ()
{
  t8_geometry_triangulated_spherical_surface *geom = new t8_geometry_triangulated_spherical_surface ();
  return (t8_geometry_c *) geom;
}

t8_geometry_c *
t8_geometry_prismed_spherical_shell_new ()
{
  t8_geometry_prismed_spherical_shell *geom = new t8_geometry_prismed_spherical_shell ();
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
t8_geometry_cubed_sphere_new ()
{
  t8_geometry_cubed_sphere *geom = new t8_geometry_cubed_sphere ();
  return (t8_geometry_c *) geom;
}

T8_EXTERN_C_END ();
