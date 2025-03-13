/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* In this test we create an occ gp_Pnt object.
 * The purpose of this test is to check whether t8code successfully links
 * against occ.
 * If t8code was not configured with --with-occ then this test
 * does nothing and is always passed.
 */

#include <t8.h>
#include <gtest/gtest.h>
#if T8_ENABLE_OCC
#include <gp_Pnt.hxx>
#endif

/* Check whether we can successfully execute VTK code */
TEST (t8_test_occ_linkage, test_gp_Pnt)
{
#if T8_ENABLE_OCC

  EXPECT_NO_THROW (gp_Pnt pnt = gp_Pnt (); pnt.SetX (1););
  t8_global_productionf ("Successfully created occ gp_Pnt object.\n");
#else
  t8_global_productionf ("This version of t8code is not compiled with occ support.\n");
#endif
}
