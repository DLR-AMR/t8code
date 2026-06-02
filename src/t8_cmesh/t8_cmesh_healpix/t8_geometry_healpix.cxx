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
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>

void
t8_geometry_healpix::t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh,
                                       [[maybe_unused]] t8_gloidx_t gtreeid, const double *ref_coords,
                                       const size_t num_coords, double *out_coords) const
{
    const t8_gloidx_t layer = gtreeid / 4;
    const t8_gloidx_t face  = gtreeid % 4;

    // Calculate the offsets based on the specific layer and face of this tree.
    const double offset_y = static_cast<double>(layer) - 1.0;
    const double offset_x = 2.0 * face + ((layer + 1) % 2);

    for (size_t i = 0; i < num_coords; ++i) {
        size_t base_idx_2d = i * 2;
        size_t base_idx_3d = i * 3; // t8code outputs 3D Cartesian coordinates (x, y, z)

        const double xi  = std::clamp(ref_coords[base_idx_2d + 0], 1e-10, 1.0 - 1e-10);
        const double eta = std::clamp(ref_coords[base_idx_2d + 1], 1e-10, 1.0 - 1e-10);

        // Transform reference coordinates into the HEALPix intermediate project space (dX, dY)
        const double dX = offset_x + (xi - eta);
        const double dY = offset_y + (xi + eta - 1.0);

        double z   = 0.0;
        double phi = 0.0;

        if (std::abs(dY) <= 1.0) {
            z   = (2.0 / 3.0) * dY;
            phi = (std::numbers::pi / 4.0) * dX;
        }
        else if (dY > 1.0) {
            z   = 1.0 - (1.0 / 3.0) * std::pow(2.0 - dY, 2.0);
            phi = (std::numbers::pi / 4.0) * (offset_x + (dX - offset_x) / (2.0 - dY));
        }
        else {
            z   = -1.0 + (1.0 / 3.0) * std::pow(2.0 + dY, 2.0);
            phi = (std::numbers::pi / 4.0) * (offset_x + (dX - offset_x) / (2.0 + dY));
        }

        z = std::clamp(z, -1.0, 1.0);

        // Convert Spherical coordinates (z, phi) to 3D Cartesian coordinates (x, y, z)
        const double theta = std::acos(z);
        const double sin_theta = std::sin(theta);

        out_coords[base_idx_3d + 0] = sin_theta * std::cos(phi);
        out_coords[base_idx_3d + 1] = sin_theta * std::sin(phi);
        out_coords[base_idx_3d + 2] = z;
    }
}
