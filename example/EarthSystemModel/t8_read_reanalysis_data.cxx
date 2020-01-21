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

#include <t8.h>
#include <sc_options.h>
#include <sc_refcount.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vec.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
#endif

/* Convert longitude and latitude coordinates to x,y,z coordinates */
void
t8_reanalysis_long_lat_to_euclid (const double longitude,
                                  const double latitude, const double R,
                                  double euclidean[3])
{
  const double        sinlong = sin (longitude);
  const double        coslong = cos (longitude);
  const double        sinlat = sin (latitude);
  const double        coslat = cos (latitude);

  euclidean[0] = R * sinlong * coslat;
  euclidean[1] = R * sinlong * sinlat;
  euclidean[2] = R * coslong;
}

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
  /* TODO: switch to T8_ALLOC */
  *pdimension_names =
    (char (*)[BUFSIZ]) malloc (number_of_dims * sizeof (*dimension_names));
  dimension_names = *pdimension_names;
  *pdimension_lengths = T8_ALLOC (size_t, number_of_dims);
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

  *pdata = T8_ALLOC (double, number_of_entries);
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

static int
t8_netcdf_open_file (const char *filename, const double radius,
                     double **pcoordinates_euclidean,
                     size_t * pnum_coordinates)
{
  int                 ncid, retval;
  int                 number_of_dims;
  int                 longitude_pos, latitude_pos;
  char                (*dimension_names)[BUFSIZ];
  size_t             *dimension_lengths;
#define NUM_DATA 3
  double             *data_in[NUM_DATA];
  double             *coordinates_euclidean;

  /* Open the file */
  t8_debugf ("Opening file %s\n", filename);
  retval = nc_open (filename, NC_NOWRITE, &ncid);
  if (retval) {
    /* Could not open the file */
    T8_NETCDF_ERROR (filename, "opening file", retval);
    return retval;
  }

  /* read the dimensions */
  retval =
    t8_netcdf_read_dimensions (filename, ncid, &number_of_dims,
                               &dimension_names, &dimension_lengths);
  if (retval) {
    /* An error occured and was printed,
     * the file is closed and we exit. */
    return retval;
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
      return retval;
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

  /* Find the position of the longitude and latitude entries */
  latitude_pos = -1;
  longitude_pos = -1;
  for (int i = 0; i < NUM_DATA; ++i) {
    if (!strcmp (dimension_names[i], "longitude")) {
      longitude_pos = i;
    }
    else if (!strcmp (dimension_names[i], "latitude")) {
      latitude_pos = i;
    }
  }
  if (latitude_pos < 0 || longitude_pos < 0) {
    t8_errorf
      ("ERROR: The fields 'longitute' and 'latitude' were not found in %s\n",
       filename);
  }
#ifdef T8_ENABLE_DEBUG
  else {
    t8_debugf ("'longitude' at pos %i, 'latitude' at pos %i\n", longitude_pos,
               latitude_pos);
  }
#endif

  /* Convert longitude and latitude to x,y,z */
  const size_t        num_long = dimension_lengths[longitude_pos];
  const size_t        num_lat = dimension_lengths[latitude_pos];
  /* Compute the number of coordinates */
  size_t              num_coordinates = *pnum_coordinates =
    num_long * num_lat;
  /* Allocate array to store all x,y,z coordinates */
  coordinates_euclidean = *pcoordinates_euclidean =
    T8_ALLOC (double, 3 * num_coordinates);

  /* Loop over all longitudes and all latitudes and compute the euclidean
   * coordinates for each point. */
  for (size_t ilong = 0; ilong < num_long; ++ilong) {
    const double        longitude = data_in[longitude_pos][ilong];
    for (size_t ilat = 0; ilat < num_lat; ++ilat) {
      const double        latitude = data_in[latitude_pos][ilat];
      double              xyz[3];
      /* Compute the current position in the euclidean array */
      const size_t        position = 3 * (num_lat * ilong + ilat);

      /* Compute euclidean coordinates of this point */
      t8_reanalysis_long_lat_to_euclid (longitude, latitude, radius, xyz);

      /* Store into the array */
      coordinates_euclidean[position] = xyz[0];
      coordinates_euclidean[position + 1] = xyz[1];
      coordinates_euclidean[position + 2] = xyz[2];
      //     t8_debugf ("%.3f %.3f %.3f - %.2f\n", xyz[0], xyz[1], xyz[2],
      //                t8_vec_norm (xyz));
    }
  }

  /* Clean-up memory */
  for (int i = 0; i < NUM_DATA; ++i) {
    T8_FREE (data_in[i]);
  }
  T8_FREE (dimension_lengths);
  free (dimension_names);

  /* Return success */
  return 0;
#undef NUM_DATA
}
#endif

/* Read msh-file and build a uniform forest on it.
 * Return 0 on success */
t8_forest_t
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
    return NULL;
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

  /* return sucess */
  return forest;
}

/* search callback function that identifies elements that contain a point.
 * This function returns true if a given point is contained in the element.
 * If the element is a leaf, we use this function to associate the element's index
 * with the point, such that we can access it later.
 * A point may be contained in multiple elements (boundaries, round-off errors,
 * or a given search tolerance (to be implemented later))
 */

typedef struct
{
  const double       *coordinates;      /* The array of coordinates of all points */
  sc_array_t         *matching_elements;        /* For each point an array of the element indices
                                                   that contain this point. (filled in the search query callback)
                                                 */
#ifdef T8_ENABLE_DEBUG
  int                 matched_elements; /* In debug mode, count how many matching elements we find. */
  int                 matched_points;   /* In debug mode, count for how many elements */
#endif
} t8_netcdf_search_user_data_t;

static int
t8_netcdf_find_mesh_elements_query (t8_forest_t forest,
                                    t8_locidx_t ltreeid,
                                    const t8_element_t *
                                    element,
                                    const int is_leaf,
                                    t8_element_array_t *
                                    leaf_elements,
                                    t8_locidx_t tree_leaf_index, void *point,
                                    size_t point_index)
{
  if (point == NULL) {
    /* The callback is called in element mode, and not in query mode.
     * We have to decide with which elements we continue the search.
     * Since we continue with every element that has matching points, we
     * do not exclude any element in this stage. Elements that do not match
     * any queries in the query stage are excluded from the search automatically */
    return 1;
  }
  double             *tree_vertices =
    t8_forest_get_tree_vertices (forest, ltreeid);
  if (t8_forest_element_point_inside
      (forest, ltreeid, element, tree_vertices, (const double *) point)) {
    /* This point is contained in this element */
    if (is_leaf) {
      /* This element is a leaf element, we add its index to the list of
       * elements that contain this point */
      t8_netcdf_search_user_data_t *user_data =
        (t8_netcdf_search_user_data_t *) t8_forest_get_user_data (forest);
#ifdef T8_ENABLE_DEBUG
      /* In debugging mode we count how many points and elements we match */
      user_data->matched_elements++;
      if (user_data->matching_elements[point_index].elem_count == 0) {
        /* This point was not found inside an element yet, we add to the counter of matched points */
        user_data->matched_points++;
      }
#endif
      /* Compute the forest local index of the element */
      t8_locidx_t         element_index =
        t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;
      /* Add this index to the array of found elements */
      *(t8_locidx_t *) sc_array_push (user_data->matching_elements +
                                      point_index) = element_index;
    }
    /* Since the point is contained in the element, we return 1 */
    return 1;
  }
  else {
    /* This point is not contained in this element, return 0 */
    return 0;
  }
}

/* Given a forest and an array of points, identify the elements that contain
 * the points.
 * The points are given as one coordinate array in the format (x_0 y_0 z_0 x_1 y_1 z_1 ... )
 */
void
t8_netcdf_find_mesh_elements (t8_forest_t forest, double *points,
                              const size_t num_points)
{
  sc_array_t         *matching_elements;
  t8_netcdf_search_user_data_t coords_and_matching_elements;
  size_t              ipoint;
  sc_array_t          queries;

  /* Allocate as many arrays as we have points to store the
   * matching elements for each point */
  matching_elements = T8_ALLOC (sc_array_t, num_points);
  /* Initialize these arrays */
  for (ipoint = 0; ipoint < num_points; ++ipoint) {
    sc_array_init (matching_elements + ipoint, sizeof (t8_locidx_t));
  }
  /* Init user data for the search routine */
  coords_and_matching_elements.coordinates = points;
  coords_and_matching_elements.matching_elements = matching_elements;
#ifdef T8_ENABLE_DEBUG
  coords_and_matching_elements.matched_elements = 0;
  coords_and_matching_elements.matched_points = 0;
#endif
  /* Set this data as the forests user pointer */
  t8_forest_set_user_data (forest, &coords_and_matching_elements);

  /* Initialize the array of points to be passed to the search function.
   * Each entry is one point, thus 3 doubles */
  sc_array_init_data (&queries, points, 3 * sizeof (double), num_points);
  t8_debugf ("Starting search with %zd points\n", num_points);
  t8_forest_search (forest, t8_netcdf_find_mesh_elements_query,
                    t8_netcdf_find_mesh_elements_query, &queries);

  t8_debugf ("Finished search. Found %i points and matched %i elements\n",
             coords_and_matching_elements.matched_points,
             coords_and_matching_elements.matched_elements);

  T8_FREE (matching_elements);
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
    t8_forest_t         forest =
      t8_reanalysis_build_forest (mesh_filename, sphere_radius, sphere_dim,
                                  comm);
    if (forest != NULL) {
      double             *coordinates_euclidean;
      size_t              num_coordinates;
      retval = t8_netcdf_open_file (netcdf_filename, sphere_radius,
                                    &coordinates_euclidean, &num_coordinates);
      if (!retval) {
        t8_netcdf_find_mesh_elements (forest, coordinates_euclidean,
                                      num_coordinates);
        /* Clean-up */
        T8_FREE (coordinates_euclidean);
      }
      t8_forest_unref (&forest);
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
