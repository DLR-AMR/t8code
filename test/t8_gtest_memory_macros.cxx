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

#include <test/t8_gtest_memory_macros.hxx>

/**
 * Package id for the testsuite. Used for attributes.
 */
static int testsuite_package_id = -1;

void
t8_testsuite_register_package_id ()
{
  /* Register a package id for the t8code testsuite */
  testsuite_package_id = sc_package_register (NULL, SC_LP_DEFAULT, "t8code_testsuite", "t8code testsuite package.");
}

int
t8_testsuite_get_package_id ()
{
  return testsuite_package_id;
}
