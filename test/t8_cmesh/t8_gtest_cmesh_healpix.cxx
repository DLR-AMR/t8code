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

/* This gtest checks whether a cmesh of a healpix is committed and whether the points have the same radius*/

#include <gtest/gtest.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_types/t8_vec.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_geometry_healpix.hxx>
#include <t8_cmesh/t8_cmesh_healpix/t8_cmesh_healpix.h>
#include <test/t8_gtest_macros.hxx>
#include <random>
/* Test whether the healpix geometry has an equal radius all across the geometry. */
TEST (T8GeometryHealpixTest, AllTreesAreOnUnitSphere)
{
  t8_geometry_healpix geom;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  std::random_device rd;
  std::mt19937 gen (rd ());
  std::uniform_real_distribution<double> dist (0.0, 1.0);

  // A set of 2D reference coordinates within the element
  std::vector<double> ref_coords = {
    0.0, 0.0,  // Corners
    1.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    0.5, 0.5,  // Center
  };
  int num_gen_coords = 0;
// in addition to the basic test cases generate random coordinates to test out implementation
#if T8_TEST_LEVEL_INT >= 2
  num_gen_coords = 10;  // Test level basic
#else
  num_gen_coords = 20;
#endif

  for (int i = 0; i < num_gen_coords; i++) {
    std::cout << i << std::endl;
    double point_x = dist (gen);
    double point_y = dist (gen);
    ref_coords.push_back (point_x);
    ref_coords.push_back (point_y);
  }

  const size_t num_coords = ref_coords.size () / 2;
  std::vector<double> out_coords (num_coords * 3, 0.0);
  // Loop through all 12 base trees of the HEALPix projection
  for (t8_gloidx_t gtreeid = 0; gtreeid < 12; ++gtreeid) {
    geom.t8_geom_evaluate (cmesh, gtreeid, ref_coords.data (), num_coords, out_coords.data ());

    for (size_t i = 0; i < num_coords; ++i) {
      std::span<const double> point (out_coords.data () + 3 * i, 3);
      double distance = t8_norm (point);
      // Expect distance to be 1.0 (Unit Sphere)
      EXPECT_NEAR (distance, 1.0, 1e-7) << "Radius check failed at Tree ID: " << gtreeid << " for point index: " << i
                                        << " (xi=" << ref_coords[i * 2] << ", eta=" << ref_coords[i * 2 + 1] << ")";
    }
  }
  t8_cmesh_destroy (&cmesh);
}
