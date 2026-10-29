/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_gtest_face_transform.cxx
 * Tests for the corner permutation of face connections, \ref t8_eclass_face_orientation_to_corner_permutation,
 * and for the face transformation of the standalone scheme, which is built on it. The default scheme is the
 * reference for the orientation convention.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <test/t8_gtest_macros.hxx>
#include <test/t8_gtest_custom_assertion.hxx>
#include <algorithm>
#include <array>
#include <map>
#include <vector>

#define T8_FACE_TRANSFORM_MAX_LEVEL 3

static const t8_eclass_t t8_face_eclasses[] = { T8_ECLASS_VERTEX, T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE };

/* The permutation is a permutation, maps corner 0 of the smaller face to the orientation,
 * and the other face maps back with its inverse. */
TEST (face_transform, corner_permutation)
{
  for (const t8_eclass_t face_eclass : t8_face_eclasses) {
    const int num_corners = t8_eclass_num_vertices[face_eclass];
    for (int orientation = 0; orientation < num_corners; orientation++) {
      for (int sign = 0; sign < 2; sign++) {
        int perm_smaller[T8_ECLASS_MAX_CORNERS_2D];
        int perm_bigger[T8_ECLASS_MAX_CORNERS_2D];
        t8_eclass_face_orientation_to_corner_permutation (face_eclass, orientation, sign, 1, perm_smaller);
        t8_eclass_face_orientation_to_corner_permutation (face_eclass, orientation, sign, 0, perm_bigger);
        EXPECT_EQ (perm_smaller[0], orientation);
        std::vector<int> sorted (perm_smaller, perm_smaller + num_corners);
        std::sort (sorted.begin (), sorted.end ());
        for (int icorner = 0; icorner < num_corners; icorner++) {
          EXPECT_EQ (sorted[icorner], icorner);
          EXPECT_EQ (perm_bigger[perm_smaller[icorner]], icorner);
        }
      }
    }
  }
}

/* The corners of an element in reference coordinates, sorted. */
static std::vector<std::array<double, 2>>
t8_test_sorted_corners (const t8_scheme *scheme, const t8_eclass_t eclass, const t8_element_t *element)
{
  std::vector<std::array<double, 2>> corners (scheme->element_get_num_corners (eclass, element));
  for (size_t icorner = 0; icorner < corners.size (); icorner++) {
    double coords[3] = { 0, 0, 0 };
    scheme->element_get_vertex_reference_coords (eclass, element, icorner, coords);
    corners[icorner] = { coords[0], coords[1] };
  }
  std::sort (corners.begin (), corners.end ());
  return corners;
}

struct face_transform_compare: public testing::TestWithParam<t8_eclass_t>
{
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    default_scheme = t8_scheme_new_default ();
    standalone_scheme = t8_scheme_new_standalone ();
  }
  void
  TearDown () override
  {
    default_scheme->unref ();
    standalone_scheme->unref ();
  }
  t8_eclass_t eclass;
  const t8_scheme *default_scheme;
  const t8_scheme *standalone_scheme;
};

/* The standalone face transformation maps every element onto the same element as the default scheme,
 * for all orientations, signs and both sides of the connection. Elements are matched by their corners,
 * since the two schemes have different space-filling curves. */
TEST_P (face_transform_compare, matches_default)
{
  t8_element_t *elem, *transformed;
  const int num_corners = t8_eclass_num_vertices[eclass];
  for (int level = 0; level <= T8_FACE_TRANSFORM_MAX_LEVEL; level++) {
    const t8_gloidx_t num_elements = default_scheme->count_leaves_from_root (eclass, level);
    ASSERT_EQ (num_elements, standalone_scheme->count_leaves_from_root (eclass, level));
    for (int orientation = 0; orientation < num_corners; orientation++) {
      for (int sign = 0; sign < 2; sign++) {
        for (int is_smaller = 0; is_smaller < 2; is_smaller++) {
          /* The image of each element of the default scheme. */
          std::map<std::vector<std::array<double, 2>>, std::vector<std::array<double, 2>>> default_images;
          default_scheme->element_new (eclass, 1, &elem);
          default_scheme->element_new (eclass, 1, &transformed);
          for (t8_gloidx_t id = 0; id < num_elements; id++) {
            default_scheme->element_set_linear_id (eclass, elem, level, id);
            /* The default quad does not set the level of a different output element for sign = 1,
             * so we transform in place, as the forest does. */
            default_scheme->element_copy (eclass, elem, transformed);
            default_scheme->element_transform_face (eclass, transformed, transformed, orientation, sign, is_smaller);
            default_images[t8_test_sorted_corners (default_scheme, eclass, elem)]
              = t8_test_sorted_corners (default_scheme, eclass, transformed);
          }
          default_scheme->element_destroy (eclass, 1, &elem);
          default_scheme->element_destroy (eclass, 1, &transformed);

          standalone_scheme->element_new (eclass, 1, &elem);
          standalone_scheme->element_new (eclass, 1, &transformed);
          for (t8_gloidx_t id = 0; id < num_elements; id++) {
            standalone_scheme->element_set_linear_id (eclass, elem, level, id);
            standalone_scheme->element_transform_face (eclass, elem, transformed, orientation, sign, is_smaller);
            EXPECT_TRUE (standalone_scheme->element_is_valid (eclass, transformed));
            EXPECT_EQ (standalone_scheme->element_get_level (eclass, transformed), level);
            const auto image = default_images.find (t8_test_sorted_corners (standalone_scheme, eclass, elem));
            ASSERT_NE (image, default_images.end ());
            EXPECT_EQ (image->second, t8_test_sorted_corners (standalone_scheme, eclass, transformed))
              << "level " << level << " id " << id << " orientation " << orientation << " sign " << sign
              << " is_smaller " << is_smaller;

            /* The transformation in place gives the same element */
            standalone_scheme->element_transform_face (eclass, elem, elem, orientation, sign, is_smaller);
            EXPECT_ELEM_EQ (standalone_scheme, eclass, elem, transformed);
          }
          standalone_scheme->element_destroy (eclass, 1, &elem);
          standalone_scheme->element_destroy (eclass, 1, &transformed);
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_face_transform, face_transform_compare,
                          testing::Values (T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE),
                          [] (const testing::TestParamInfo<t8_eclass_t> &info) {
                            return std::string (t8_eclass_to_string[info.param]);
                          });
