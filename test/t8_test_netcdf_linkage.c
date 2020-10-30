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

#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
#endif

static void
t8_test_netcdf_linkage ()
{
#if T8_WITH_NETCDF

/* Create an in-memory netcdf file. This file will not be stored on
 * the disk. */
  int                 ncid;
  int                 nc_error = nc_create ("FileName", NC_DISKLESS, &ncid);

  /* Check for error */
  SC_CHECK_ABORTF (nc_error == NC_NOERR,
                   "netcdf error when creating in memory " "file: %s\n",
                   nc_strerror (nc_error));

  /* Close the file */
  nc_error = nc_close (ncid);
  /* Check for error */
  SC_CHECK_ABORTF (nc_error == NC_NOERR,
                   "netcdf error when closing in memory " "file: %s\n",
                   nc_strerror (nc_error));

  t8_global_productionf
    ("Successfully created and closed in memory netcdf file.\n");
#else
  t8_global_productionf
    ("This version of t8code is not compiled with netcdf support.\n");
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  /* Check netcdf linkage */
  t8_test_netcdf_linkage ();

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
