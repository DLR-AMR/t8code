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

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_readmshfile.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
#endif

/* TODO: Use T8_ALLOC instead of malloc */

#if T8_WITH_NETCDF
/* Print an error message using the netcdf errorcode */
#define T8_NETCDF_ERROR(filename, description, errorcode) \
  t8_errorf("Error in file %s - %s - %s\n", filename, description, nc_strerror(errorcode))

/* Close an opened netcdf file */
static int
t8_netcdf_close_file (const char *filename, int ncid)
{
  int                 retval;

  /* Close the file, freeing all resources. */
  t8_debugf ("Closing file %s\n", filename);
  retval = nc_close (ncid);
  if (retval) {
    /* Could not close the file */
    T8_NETCDF_ERROR (filename, "closing file", retval);
  }
  return retval;
}

/**
 * From an opened netcdf file read the names and lengths of the stored dimensions.
 * \param[in] filename            The filename of the opened netcdf file (used only for error output)
 * \param[in] ncid                The file id
 * \param[out] pnumber_of_dims    On output the number of dimensions in the file
 * \param[out] pdimension_names   On output allocated to \a *pnumber_of_dims entries and stores the names of the dimensions
 * \param[out] pdimension_lengths On output allocated to \a *pnumber_of_dims entries and stores the lengths of the dimensions
 * \return                        netcdf error value. 0 on success.
 */
static int
t8_netcdf_read_dimensions (const char *filename, const int ncid,
                           int *pnumber_of_dims,
                           char (**pdimension_names)[BUFSIZ],
                           size_t ** pdimension_lengths)
{
  int                 dimension_id;
  int                 retval;
  int                 number_of_dims;
  char                (*dimension_names)[BUFSIZ];
  size_t             *dimension_lengths;

  /* read the number of dimension ids */
  retval = nc_inq_ndims (ncid, pnumber_of_dims);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading number of dimensions", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }
  number_of_dims = *pnumber_of_dims;
  /* Allocate dimension_names and length arrays */
  *pdimension_names =
    (char (*)[BUFSIZ]) malloc (number_of_dims * sizeof (*dimension_names));
  dimension_names = *pdimension_names;
  *pdimension_lengths =
    (size_t *) malloc (number_of_dims * sizeof (*pdimension_lengths));
  dimension_lengths = *pdimension_lengths;
  if (dimension_names == NULL || dimension_lengths == NULL) {
    t8_global_errorf ("Could not allocate memory for %i dimension names\n",
                      number_of_dims);
    t8_netcdf_close_file (filename, ncid);
    free (dimension_names);
    free (dimension_lengths);
    return retval;
  }
  t8_debugf ("Reading %i dimensions...\n", number_of_dims);

  /* Read the names and length of the dimensions */
  for (dimension_id = 0; dimension_id < number_of_dims; ++dimension_id) {
    retval =
      nc_inq_dim (ncid, dimension_id, dimension_names[dimension_id],
                  dimension_lengths + dimension_id);
    if (retval) {
      T8_NETCDF_ERROR (filename, "reading dimension names and lengths",
                       retval);
      t8_netcdf_close_file (filename, ncid);
      return retval;
    }
    t8_debugf ("Read dimension [%s] of length %zu\n",
               dimension_names[dimension_id],
               dimension_lengths[dimension_id]);
  }
  /* return success */
  return 0;
}

static int
t8_netcdf_read_double_data (const char *filename, const int ncid,
                            const char *varname,
                            const int expected_number_of_dims,
                            const int number_of_entries, double **pdata)
{
  int                 varid;
  int                 dimension_id;
  int                 retval;
  double             *data;

  /* Get the varid of the longitude data variable, based on its name. */
  t8_debugf ("Reading data info for '%s'\n", varname);
  retval = nc_inq_varid (ncid, varname, &varid);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading data info", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

#ifdef T8_ENABLE_DEBUG
  /* Ensure that the variable has the expected number of dimensions */
  {
    int                 ndims;
    retval = nc_inq_varndims (ncid, varid, &ndims);
    if (retval) {
      T8_NETCDF_ERROR (filename, "reading number of dimensions", retval);
    }
    if (ndims != expected_number_of_dims) {
      t8_global_errorf
        ("Error: '%s' variable has more than %i dimension\n",
         varname, expected_number_of_dims);
      t8_netcdf_close_file (filename, ncid);
      return retval;
    }
    else {
      t8_debugf ("'%s' has exactly 1 dimension as expected\n", varname);
    }
  }
#endif

  /* Read the number of entries */
  retval = nc_inq_vardimid (ncid, varid, &dimension_id);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading dimension id", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  *pdata = (double *) malloc (number_of_entries * sizeof (**pdata));
  data = *pdata;
  if (data == NULL) {
    t8_global_errorf ("Could not allocate memory for %i data items\n",
                      number_of_entries);
    t8_netcdf_close_file (filename, ncid);
  }
  t8_debugf ("'%s' has %i entries\n", varname, number_of_entries);

  /* Read the longitude data. */
  t8_debugf ("Reading '%s' data\n", varname);
  retval = nc_get_var_double (ncid, varid, data);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading data", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }
  /* return success */
  return 0;
}

static void
t8_netcdf_open_file (const char *filename)
{
  int                 ncid, retval;
  int                 number_of_dims;
  char                (*dimension_names)[BUFSIZ];
  size_t             *dimension_lengths;
#define NUM_DATA 3
  double             *data_in[NUM_DATA];

  /* Open the file */
  t8_debugf ("Opening file %s\n", filename);
  retval = nc_open (filename, NC_NOWRITE, &ncid);
  if (retval) {
    /* Could not open the file */
    T8_NETCDF_ERROR (filename, "opening file", retval);
    return;
  }

  /* read the dimensions */
  retval =
    t8_netcdf_read_dimensions (filename, ncid, &number_of_dims,
                               &dimension_names, &dimension_lengths);
  if (retval) {
    /* An error occured and was printed,
     * the file is closed and we exit. */
    return;
  }

  if (number_of_dims < NUM_DATA) {
    t8_global_errorf ("Expected at least 3 dimensions. Only %i found.\n",
                      number_of_dims);
    t8_netcdf_close_file (filename, ncid);
  }
  for (int i = 0; i < NUM_DATA; ++i) {
    retval =
      t8_netcdf_read_double_data (filename, ncid, dimension_names[i], 1,
                                  dimension_lengths[i], &data_in[i]);
    if (retval) {
      /* An error occured and was printed,
       * the file is closed and we exit. */
      return;
    }
  }

  /* Close the opened file */
  t8_netcdf_close_file (filename, ncid);

#ifdef T8_ENABLE_DEBUG
  /* Output the read data */
  {
    for (int i = 0; i < NUM_DATA; ++i) {
      size_t              j;
      char                output[BUFSIZ] = "";
      char                number[20];
      for (j = 0; j < dimension_lengths[i]; ++j) {
        snprintf (number, 20, " %.2f", data_in[i][j]);
        if (strlen (output) < BUFSIZ - 21) {
          strcat (output, number);
        }
      }
      t8_debugf ("Read data from '%s' with %zd entries:\n",
                 dimension_names[i], dimension_lengths[i]);
      t8_debugf ("%s\n", output);
    }
  }
#endif

  /* Clean-up memory */
  for (int i = 0; i < NUM_DATA; ++i) {
    free (data_in[i]);
  }
  free (dimension_lengths);
  free (dimension_names);
#undef NUM_DATA
}
#endif

/* Read msh-file and build a uniform forest on it.
 * Return 0 on success */
int
t8_reanalysis_build_forest (const char *mesh_filename, double radius,
                            int dimension, sc_MPI_Comm comm)
{
  const int           level = 0;
  const int           do_ghosts = 1;
  /* read the coarse mesh from the .msh file */
  t8_cmesh_t          cmesh =
    t8_cmesh_from_msh_file (mesh_filename, 0, comm, dimension, 0);
  if (cmesh == NULL) {
    /* cmesh could not be built */
    t8_global_errorf ("Error when openening file %s\n", mesh_filename);
    return 1;
  }
  /* build a uniform forest from the coarse mesh */
  t8_forest_t         forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (),
                           level, do_ghosts, comm);

#ifdef T8_ENABLE_DEBUG
  /* in debug mode, write the forest to vtk */
  {
    char                output_file_prefix[BUFSIZ];
    snprintf (output_file_prefix, BUFSIZ, "forest_uniform_l%i_%s",
              level, mesh_filename);
    t8_forest_write_vtk (forest, output_file_prefix);
  }
#endif

  t8_forest_unref (&forest);

  /* return sucess */
  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  const char         *netcdf_filename = NULL;
  const char         *mesh_filename = NULL;
  int                 parsed, helpme;
  int                 sphere_dim;
  double              sphere_radius;
  const sc_MPI_Comm   comm = sc_MPI_COMM_WORLD; /* The mpi communicator used throughout */

  /* help message, prints when called with '-h' option */
  snprintf (help, BUFSIZ, "This program reads data from a netcdf file.\n");

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize libsc */
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
#ifdef T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_PRODUCTION);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  /* Add command line arguments */
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_string (opt, 'n', "netcdffile", &netcdf_filename, NULL,
                         "The netcdf-file that should be read.");
  sc_options_add_string (opt, 'f', "meshfile", &mesh_filename, NULL,
                         "The msh-file of a sphere that should be read (without the '.msh').");
  sc_options_add_double (opt, 'r', "radius", &sphere_radius, 1.0,
                         "The radius of the sphere in the msh file. Default = 1");
  sc_options_add_int (opt, 'd', "dim", &sphere_dim, 2,
                      "The dimension of the mesh. Default = 2");

  /* Parse the command line arguments from the input */
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
#if T8_WITH_NETCDF
  else if (parsed >= 0 && netcdf_filename != NULL && mesh_filename != NULL) {
    int                 retval;
    retval =
      t8_reanalysis_build_forest (mesh_filename, sphere_radius, sphere_dim,
                                  comm);
    if (!retval) {
      t8_netcdf_open_file (netcdf_filename);
    }
  }
  else {
    /* Error when parsing the arguments */
    /* wrong usage */
    t8_global_essentialf ("\n\tERROR: Wrong usage.\n\n");
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
#else
  /* If t8code is not linked against netcdf, we cannot execute the code */
  /* Print help message and exit */
  t8_global_essentialf ("%s\n", help);
  t8_global_essentialf ("t8code is not linked against netcdf.\n");
  t8_global_essentialf
    ("To run this example configure t8code with '--with-netcdf'.\n");
#endif

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
