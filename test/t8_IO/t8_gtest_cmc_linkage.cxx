/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/* The purpose of this test is to check whether t8code successfully links
 * against cmc.
 * If t8code was not compiled with -DT8CODE_ENABLE_CMC=ON then this test
 * does nothing and is always passed.
 */
#include <gtest/gtest.h>
#include <t8.h>
#if T8_ENABLE_CMC
#include <cmc.hxx>
#include <utilities/cmc_log_functions.hxx>
#endif

/* Test correct macro dependencies.
 * Will throw a compile time error if T8_VTK_VERSION_USED
 * is defined but T8_ENABLE_VTK is not. */
#ifndef T8_ENABLE_CMC
#ifdef T8_CMC_VERSION_USED
#error Configuration error: T8_CMC_VERSION_USED is defined despite \
 T8_ENABLE_CMC not being defined.
#endif
#endif

/* Check whether T8_VTK_VERSION_USED equals VTK_MAJOR_VERSION.VTK_MINOR_VERSION */
TEST (t8_gtest_cmc_linkage, t8_test_cmc_linkage)
{
#if T8_ENABLE_CMC
  t8_debugf ("Using cmc version %s.\n", CMC_VERSION);
  cmc::CmcInitialize(cmc::kMinimumInitialization);
  cmc::cmc_global_msg("Global message via cmc.");
  cmc::CmcFinalize();
#else
  t8_debugf ("This version of t8code is not compiled with cmc support.\n");
#endif
}
