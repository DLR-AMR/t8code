/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

/** \file t8_gtest_geometry_macros.hxx
 * This file contains macros and base classes for geometry tests.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>

#ifdef T8_ENABLE_LESS_TESTS
#define T8_NUM_SAMPLE_POINTS 1000
#else
#define T8_NUM_SAMPLE_POINTS 10000
#endif

class geometry_test: public testing::TestWithParam<t8_eclass_t> {
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
  }
  void
  TearDown () override
  {
  }
  t8_eclass_t eclass;
};
