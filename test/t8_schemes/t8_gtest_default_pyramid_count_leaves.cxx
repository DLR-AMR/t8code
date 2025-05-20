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

#include <gtest/gtest.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

TEST (class_element_leaves, test_default_pyramid_count_leaves_root)
{

  const t8_scheme *scheme;
  scheme = t8_scheme_new_default ();
  t8_eclass_t eclass = T8_ECLASS_PYRAMID;

  const int maxlevel = scheme->get_maxlevel (eclass);
  t8_gloidx_t compare_value = 1;
  t8_gloidx_t sum1 = 1;
  t8_gloidx_t sum2 = 1;

  for (int level = 0; level <= maxlevel; ++level) {
    const t8_gloidx_t leaf_count = scheme->count_leaves_from_root (eclass, level);
    ASSERT_EQ (leaf_count, compare_value)
      << "Incorrect leaf count " << leaf_count << " at eclass " << t8_eclass_to_string[eclass] << " and level " << level
      << " (expecting " << compare_value << ")";
    //Calculate the copare value for the default pyramid
    sum1 *= 8;
    sum2 *= 6;
    compare_value = 2 * sum1 - sum2;
  }
  scheme->unref ();
}
