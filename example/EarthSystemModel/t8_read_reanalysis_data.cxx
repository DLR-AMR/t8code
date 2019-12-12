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

/* TOOD: Need T8_WITH_NETCDF macro? */
/* TODO: Use T8_ALLOC instead of malloc */

/* Print an error message using the netcdf errorcode */
#define T8_NETCDF_ERROR(filename, description, errorcode) \
  t8_debugf("Error in file %s - %s - %s\n", filename, description, nc_strerror(errorcode))

/* Close an opened netcdf file */
static void
closeFile (const char *filename, int ncid)
{
  int                 retval;

  /* Close the file, freeing all resources. */
  t8_debugf ("Closing file %s\n", filename);
  retval = nc_close (ncid);
  if (retval) {
    /* Could not close the file */
    T8_NETCDF_ERROR (filename, "closing file", retval);
  }
}

static void
openFile ()
{
  /* This currently fails since the file does not exist yet.
   * We currently use this call to calibrate the linker settings
   * for netcdf. */
  int                 ncid, retval, varid;
  const char          filename[BUFSIZ] = "./July2019_12_SurfaceUVSnow.nc";
#define NUM_LONGITUDE 480
  int                 number_of_entries = 0;
  int                 number_of_dims;
  char                (*dimension_names)[BUFSIZ];
  size_t             *dimension_length;
  int                 dimension_id;
  double             *data_in;  /* TODO: Read number of entries from file. Currently hardcoded */

  /* Open the file */
  t8_debugf ("Opening file %s\n", filename);
  retval = nc_open (filename, NC_NOWRITE, &ncid);
  if (retval) {
    /* Could not open the file */
    T8_NETCDF_ERROR (filename, "opening file", retval);
    return;
  }

  /* read the number of dimension ids */
  retval = nc_inq_ndims (ncid, &number_of_dims);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading number of dimensions", retval);
    closeFile (filename, ncid);
    return;
  }
  /* Allocate dimension_names and lenght arrays */
  dimension_names =
    (char (*)[BUFSIZ]) malloc (number_of_dims * sizeof (*dimension_names));
  dimension_length =
    (size_t *) malloc (number_of_dims * sizeof (*dimension_length));
  if (dimension_names == NULL || dimension_length == NULL) {
    t8_global_errorf ("Could not allocate memory for %i dimension names\n",
                      number_of_dims);
    closeFile (filename, ncid);
    free (dimension_names);
    free (dimension_length);
    return;
  }
  t8_debugf ("Reading %i dimensions...\n", number_of_dims);

  /* Read the names and length of the dimensions */
  for (dimension_id = 0; dimension_id < number_of_dims; ++dimension_id) {
    retval =
      nc_inq_dim (ncid, dimension_id, dimension_names[dimension_id],
                  dimension_length + dimension_id);
    if (retval) {
      T8_NETCDF_ERROR (filename, "reading dimension names and lengths",
                       retval);
      closeFile (filename, ncid);
      return;
    }
    t8_debugf ("Read dimension [%s] of length %zu\n",
               dimension_names[dimension_id], dimension_length[dimension_id]);
  }

  /* Get the varid of the longitude data variable, based on its name. */
  t8_debugf ("Reading longitude info\n");
  retval = nc_inq_varid (ncid, "longitude", &varid);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading longitude info", retval);
    closeFile (filename, ncid);
    return;
  }

#ifdef T8_ENABLE_DEBUG
  /* Ensure that the variable has exaclty one dimension */
  {
    int                 ndims;
    retval = nc_inq_varndims (ncid, varid, &ndims);
    if (retval) {
      T8_NETCDF_ERROR (filename, "reading number of dimensions", retval);
    }
    if (ndims != 1) {
      t8_global_errorf
        ("Error: longitude variable has more than 1 dimension\n");
      closeFile (filename, ncid);
      return;
    }
    else {
      t8_debugf ("longitude has exactly 1 dimension as expected\n");
    }
  }
#endif

  /* Read the number of entries */
  retval = nc_inq_vardimid (ncid, varid, &dimension_id);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading dimension id", retval);
    closeFile (filename, ncid);
    return;
  }

  number_of_entries = dimension_length[dimension_id];
  data_in = (double *) malloc (number_of_entries * sizeof (*data_in));
  if (data_in == NULL) {
    t8_global_errorf ("Could not allocate memory for %i data items\n",
                      number_of_entries);
    closeFile (filename, ncid);
  }
  t8_debugf ("longitude has %i entries\n", number_of_entries);

  /* Read the longitude data. */
  t8_debugf ("Reading longitude data\n");
  retval = nc_get_var_double (ncid, varid, data_in);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading longitude data", retval);
    closeFile (filename, ncid);
    return;
  }

  /* Close the opened file */
  closeFile (filename, ncid);

#ifdef T8_ENABLE_DEBUG
  /* Output the read data */
  {
    int                 i;
    char                output[BUFSIZ] = "";
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

  /* Clean-up memory */
  free (data_in);
  free (dimension_length);
  free (dimension_names);
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
