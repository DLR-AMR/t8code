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

#include <gtest/gtest.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_vector/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_cmesh_healpix.h>
#include <test/t8_gtest_macros.hxx>
#include <numbers>

struct t8_location_result
{
  double z;
  double phi;
  double sintheta;
  bool have_sintheta;
};

constexpr t8_location_result
xyf2loc (const double x, const double y, const int face, const int JRLL[], const int JPLL[])
{
  /*
   A helper fuction translated from the HEALPix julia code base.
   */
  double z = 0.0;
  double nr = 0.0;
  double sintheta = 0.0;
  bool have_sintheta = false;

  const double jr = static_cast<double> (JRLL[face]) - x - y;

  if (jr < 1.0) {
    nr = jr;
    const double tmp = (nr * nr) / 3.0;
    z = 1.0 - tmp;
    if (z > 0.99) {
      sintheta = std::sqrt (tmp * (2.0 - tmp));
      have_sintheta = true;
    }
  }
  else if (jr > 3.0) {
    nr = 4.0 - jr;
    const double tmp = (nr * nr) / 3.0;
    z = tmp - 1.0;
    if (z < -0.99) {
      sintheta = std::sqrt (tmp * (2.0 - tmp));
      have_sintheta = true;
    }
  }
  else {
    nr = 1.0;
    z = (2.0 - jr) * 2.0 / 3.0;
  }

  double tmp = static_cast<double> (JPLL[face]) * nr + x - y;
  if (tmp < 0.0)
    tmp += 8.0;
  if (tmp >= 8.0)
    tmp -= 8.0;

  const double phi = (std::abs (nr) < 1e-15) ? 0.0 : ((std::numbers::pi / 4.0 * tmp) / nr);

  return { z, phi, sintheta, have_sintheta };
}

TEST (t8_healpix_geometry_test, test_batch_against_xyf2loc_reference)
{
  const int JRLL[12] = { 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
  const int JPLL[12] = { 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };

  std::array<double, 7> test_samples = { 0.0, 0.001, 0.25, 0.5, 0.75, 0.999, 1.0 };
  t8_geometry_healpix geom;

  // Build a batch containing all combinations of test_samples for a single evaluation call
  std::vector<double> ref_coords;
  for (double xi : test_samples) {
    for (double eta : test_samples) {
      ref_coords.push_back (xi);
      ref_coords.push_back (eta);
    }
  }
  const size_t num_coords = test_samples.size () * test_samples.size ();

  // Loop through all t8code layers and faces
  for (int layer = 0; layer < 3; ++layer) {
    for (int face = 0; face < 4; ++face) {
      // t8code Tree ID (0..11)
      const t8_gloidx_t gtreeid = layer * 4 + face;

      // HEALPix base face index corresponding to t8code's layer
      // layer 0 (South) -> HEALPix faces 8..11
      // layer 1 (Equator) -> HEALPix faces 4..7
      // layer 2 (North) -> HEALPix faces 0..3
      const int healpix_face = (2 - layer) * 4 + face;

      // Output buffers
      std::vector<double> t8_out (num_coords * 3);

      // Call the batch evaluation method provided
      geom.t8_geom_evaluate (nullptr, gtreeid, ref_coords.data (), num_coords, t8_out.data ());

      // Validate each point in the batch against xyf2loc reference
      for (size_t i = 0; i < num_coords; ++i) {
        const double xi = ref_coords[i * 2 + 0];
        const double eta = ref_coords[i * 2 + 1];

        // Clamped values used inside t8_geom_evaluate
        const double xi_clamped = std::clamp (xi, 1e-10, 1.0 - 1e-10);
        const double eta_clamped = std::clamp (eta, 1e-10, 1.0 - 1e-10);

        // Pass the 0..11 global face index to xyf2loc
        const t8_location_result ref_loc = xyf2loc (xi_clamped, eta_clamped, healpix_face, JRLL, JPLL);

        const double ref_sin_theta
          = ref_loc.have_sintheta ? ref_loc.sintheta : std::sqrt (std::max (0.0, 1.0 - ref_loc.z * ref_loc.z));
        const double ref_x = ref_sin_theta * std::cos (ref_loc.phi);
        const double ref_y = ref_sin_theta * std::sin (ref_loc.phi);
        const double ref_z = ref_loc.z;

        const size_t idx3d = i * 3;

        // Allow a slightly relaxed tolerance near the poles where sin(theta) ~ 0
        const double adaptive_tol = (std::abs (ref_z) > 0.99) ? T8_PRECISION_SQRT_EPS : (4.0 * T8_PRECISION_EPS);
        EXPECT_NEAR (t8_out[idx3d + 0], ref_x, adaptive_tol)
          << "Mismatch in X | Tree ID: " << gtreeid << " (Layer: " << layer << ", Face: " << face << ")"
          << ", (xi, eta): (" << xi << ", " << eta << ")";
        EXPECT_NEAR (t8_out[idx3d + 1], ref_y, adaptive_tol)
          << "Mismatch in Y | Tree ID: " << gtreeid << " (Layer: " << layer << ", Face: " << face << ")"
          << ", (xi, eta): (" << xi << ", " << eta << ")";
        EXPECT_NEAR (t8_out[idx3d + 2], ref_z, adaptive_tol)
          << "Mismatch in Z | Tree ID: " << gtreeid << " (Layer: " << layer << ", Face: " << face << ")"
          << ", (xi, eta): (" << xi << ", " << eta << ")";
      }
    }
  }
}
