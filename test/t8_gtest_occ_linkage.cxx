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
 * If t8code was not configured with -DT8CODE_ENABLE_OCC=ON then this test
 * does nothing and is always passed.
 */

#include <t8.h>
#include <gtest/gtest.h>
#if T8_ENABLE_OCC
#include <gp_Pnt.hxx>
#include <Standard_Version.hxx>
#endif

/* Test correct macro dependencies.
 * Will throw a compile time error if T8_OCC_VERSION_USED
 * is defined but T8_ENABLE_OCC is not. */
#ifndef T8_ENABLE_OCC
#ifdef T8_OCC_VERSION_USED
#error Configuration error: T8_OCC_VERSION_USED is defined despite \
 T8_ENABLE_OCC not being defined.
#endif
#endif

/* Check whether T8_OCC_VERSION_USED equals OCC_VERSION_STRING */
TEST (t8_gtest_occ_linkage, t8_test_occ_version_number)
{
#if T8_ENABLE_OCC
  /* Compare the macro generated at configuration time with the linked runtime string */
  ASSERT_FALSE (strcmp (T8_OCC_VERSION_USED, OCC_VERSION_STRING))
    << "Linked OCC version (" << OCC_VERSION_STRING << ") does not equal the version t8code was configured with ("
    << T8_OCC_VERSION_USED << ").\n";

  if (!strcmp (T8_OCC_VERSION_USED, OCC_VERSION_STRING)) {
    t8_debugf ("Using OCC version %s.\n", OCC_VERSION_STRING);
  }
#endif
}

/* Check whether we can successfully execute OCC code */
TEST (t8_gtest_occ_linkage, test_gp_Pnt)
{
#if T8_ENABLE_OCC
  EXPECT_NO_THROW ({
    gp_Pnt pnt = gp_Pnt ();
    pnt.SetX (1);
  });
  t8_debugf ("Successfully created occ gp_Pnt object.\n");
#else
  t8_debugf ("This version of t8code is not compiled with occ support.\n");
#endif
}
