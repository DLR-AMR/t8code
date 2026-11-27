/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

#include <t8.h>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>
#include <numbers>

/**
 * Evaluates a point on a sphere surface geometry from mesh layer and face parameters.
 *
 * Maps 3D parametric reference coordinates (\p xi, \p eta) of a specific mesh element,
 * identified by its \p layer and \p face indices, to 3D Cartesian coordinates on a unit sphere.
 *
 * @param[in]  layer     The layer index of the mesh element.
 * @param[in]  face      The face index of the mesh element within the layer.
 * @param[in]  xi        The reference coordinate in the first parametric direction (typically [0, 1]).
 * @param[in]  eta       The reference coordinate in the second parametric direction (typically [0, 1]).
 * @param[in]  radial_pos Normalized radial position in [0.0, 1.0] (0.0 -> inner_radius, 1.0 -> outer_radius).
 * @param[out] out_coord Pointer to an array of size 3 where the resulting (x, y, z)
 *                       Cartesian coordinates will be stored.
 * @param[in]  inner_radius Inner radius of the shell (clamped to [0.0, outer_radius]).
 * @param[in]  outer_radius Outer radius of the shell.
 */
void
t8_eval_geom_point (t8_gloidx_t layer, t8_gloidx_t face, double xi, double eta, double radial_pos, double *out_coord,
                    double inner_radius, double outer_radius)
{
  const double safe_outer = std::max (0.0, outer_radius);
  const double safe_inner = std::clamp (inner_radius, 0.0, safe_outer);
  const double clamped_pos = std::clamp (radial_pos, 0.0, 1.0);

  // Interpolate between inner_radius and outer_radius; as safe_inner -> safe_outer, shell becomes thin
  const double r = safe_inner + clamped_pos * (safe_outer - safe_inner);

  const double offset_y = 1.0 - static_cast<double> (layer);
  const double offset_x = 2.0 * face + ((layer + 1) % 2);

  const double dX = offset_x + (xi - eta);
  const double dY = offset_y + (xi + eta - 1.0);

  double z = 0.0;
  double phi = 0.0;

  if (std::abs (dY) <= 1.0) {
    z = (2.0 / 3.0) * dY;
    phi = (std::numbers::pi / 4.0) * dX;
  }
  else if (dY > 1.0) {
    const double denominator = 2.0 - dY;
    z = 1.0 - (1.0 / 3.0) * std::pow (denominator, 2.0);

    if (denominator > 0.0) {
      phi = (std::numbers::pi / 4.0) * (offset_x + (dX - offset_x) / denominator);
    }
    else {
      // At the north pole, phi is arbitrary.
      phi = (std::numbers::pi / 4.0) * offset_x;
    }
  }
  else {
    const double denominator = 2.0 + dY;
    z = -1.0 + (1.0 / 3.0) * std::pow (denominator, 2.0);

    if (denominator > 0.0) {
      phi = (std::numbers::pi / 4.0) * (offset_x + (dX - offset_x) / denominator);
    }
    else {
      // At the south pole, phi is arbitrary.
      phi = (std::numbers::pi / 4.0) * offset_x;
    }
  }

  z = std::clamp (z, -1.0, 1.0);

  // Convert unit-sphere spherical coordinates (z, phi) to Cartesian coordinates at radius r
  const double theta = std::acos (z);
  const double sin_theta = std::sin (theta);

  out_coord[0] = (sin_theta * std::cos (phi)) * r;
  out_coord[1] = (sin_theta * std::sin (phi)) * r;
  out_coord[2] = z * r;
}

void
t8_geometry_healpix::t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                       const size_t num_coords, double *out_coords) const
{
  // This geometry is inspired by HEALPix
  const t8_gloidx_t layer = gtreeid / 4;
  const t8_gloidx_t face = gtreeid % 4;

  for (size_t i = 0; i < num_coords; ++i) {
    const size_t base_idx_ref = i * 3;  // 3D reference coordinates (xi, eta, zeta)
    const size_t base_idx_3d = i * 3;   // 3D Cartesian coordinates (x, y, z)

    const double xi = ref_coords[base_idx_ref + 0];
    const double eta = ref_coords[base_idx_ref + 1];
    const double zeta = ref_coords[base_idx_ref + 2];  // Radial parameter in [0, 1]

    t8_eval_geom_point (layer, face, xi, eta, zeta, out_coords + base_idx_3d, sphere_inner_radius, sphere_outer_radius);
  }
}
t8_geometry_type_t
t8_geometry_healpix::t8_geom_get_type () const
{
  return T8_GEOMETRY_TYPE_ANALYTIC;  // Or the specific t8_geometry_type_t enum value for HEALPix
}
