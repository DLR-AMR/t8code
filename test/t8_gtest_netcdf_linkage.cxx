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

/* In this test we create a netcdf in memory file and close it.
 * The purpose of this test is to check whether t8code successfully links
 * against netcdf.
 * If t8code was not configured with --with-netcdf then this test
 * does nothing and is always passed.
 */

#include <gtest/gtest.h>
#include <t8.h>
#if T8_ENABLE_NETCDF
#include <netcdf.h>
#endif

TEST (t8_gtest_netcdf_linkage, test_linking_with_netcdf)
{
#if T8_ENABLE_NETCDF

  /* Create an in-memory netcdf file. This file will not be stored on
 * the disk. */
  int ncid;
  int nc_error = nc_create ("FileName", NC_DISKLESS, &ncid);

  /* Check for error */
  ASSERT_EQ (nc_error, NC_NOERR) << "netcdf error when creating in memory file: " << nc_strerror (nc_error);

  /* Close the file */
  nc_error = nc_close (ncid);
  /* Check for error */
  ASSERT_EQ (nc_error, NC_NOERR) << "netcdf error when closing in memory file: " << nc_strerror (nc_error);

  t8_debugf ("Successfully created and closed in memory netcdf file.\n");
#else
  t8_debugf ("This version of t8code is not compiled with netcdf support.\n");
#endif
}
