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

#include <netcdf.h>
#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>

#define T8_NETCDF_ERROR(filename, description, errorcode) \
  t8_debugf("Error in file %s - %s - %s\n", filename, description, nc_strerror(errorcode))

static void
openFile ()
{
  /* This currently fails since the file does not exist yet.
   * We currently use this call to calibrate the linker settings
   * for netcdf. */
  int                 ncid, retval, varid;
  const char          filename[BUFSIZ] = "./July2019_12_SurfaceUVSnow.nc";
#define NUM_LONGITUDE 480
  double              data_in[NUM_LONGITUDE];   /* TODO: Read number of entries from file. Currently hardcoded */

  /* Open the file */
  t8_debugf ("Opening file %s\n", filename);
  retval = nc_open (filename, NC_NOWRITE, &ncid);
  if (retval) {
    /* Could not open the file */
    T8_NETCDF_ERROR (filename, "opening file", retval);
    return;
  }

  /* Get the varid of the longitude data variable, based on its name. */
  t8_debugf ("Reading longitude info\n");
  retval = nc_inq_varid (ncid, "longitude", &varid);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading longitude info", retval);
    goto closeAndExit;
  }

  /* Read the longitude data. */
  t8_debugf ("Reading longitude data\n");
  retval = nc_get_var_double (ncid, varid, data_in);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading longitude data", retval);
    goto closeAndExit;
  }

#ifdef T8_ENABLE_DEBUG
  /* Output the read data */
  {
    int                 i;
    char                output[BUFSIZ];
    char                number[20];
    for (i = 0; i < NUM_LONGITUDE; ++i) {
      /* TODO: Use more sophisticated version to append to the string
       *       Using the same buffer as in and output in snprintf causes undefined behaviour */
      snprintf (number, 20, " %.2f", data_in[i]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
    }
    t8_debugf ("%s\n", output);
  }
#endif

closeAndExit:
  /* Close the file, freeing all resources. */
  t8_debugf ("Closing file %s\n", filename);
  retval = nc_close (ncid);
  if (retval) {
    /* Could not close the file */
    T8_NETCDF_ERROR (filename, "closing file", retval);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  openFile ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
