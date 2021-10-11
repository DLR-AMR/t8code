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
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest_vtk.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vec.h>
#include <example/common/t8_example_common.h>
#if !T8_WITH_NETCDF
#else
#include <netcdf.h>

/* TODO: - Document properly
 *       - move some netcdf routines to internal netcdf header
 *       - Point search needs a range in which to look, since the geometry of
 *         the forest does not exactly match the input geometry -> tolerance of point search as parameter
 *       - Do not store the point coordinates explicitely but pass longitude and
 *         latitude to the search and compute the point coordinates when needed.
 *       - u10 and v10 veocities are tangential to longitude and latitude, we need to transform these
 *         into cartesian bases if we want to use our advection solver
 *       - use u10 and v10 as vector field for advection solver
 *       - distribute the mesh (before search)
 *       - for each element store a list of the points that match it, refine the element and
 *         search again for the children and store the new points (use forest_replace)
 *       - Refine the mesh where the gradient of u10/v10 is large
 *       - Use search to find the point with the minimal distance to an element?
 *       - Reduce File i/o to single process where possible
 *       - Add tolerance to search_data
 *       - Stop refining if we copied the data in the step before.
 *       - Refine when v10/u10 is large
 *       - seperate data loading and data to mesh
 */

/* Convert longitude and latitude coordinates to x,y,z coordinates */
void
t8_reanalysis_long_lat_to_euclid (const double longitude,
                                  const double latitude, const double R,
                                  double euclidean[3])
{
  T8_ASSERT (-90 <= latitude && latitude <= 90);
  T8_ASSERT (0 <= longitude && longitude <= 360);
  const double        sinlong = sin (longitude / 180 * M_PI);
  const double        coslong = cos (longitude / 180 * M_PI);
  const double        sinlat = sin (latitude / 180 * M_PI);
  const double        coslat = cos (latitude / 180 * M_PI);

  euclidean[0] = R * coslat * coslong;
  euclidean[1] = R * coslat * sinlong;
  euclidean[2] = R * sinlat;
}

/* TODO: Use T8_ALLOC instead of malloc */

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
    T8_FREE (dimension_lengths);
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
      free (dimension_names);
      T8_FREE (dimension_lengths);
      return retval;
    }
    t8_debugf ("Read dimension [%s] of length %zu\n",
               dimension_names[dimension_id],
               dimension_lengths[dimension_id]);
  }
  /* return success */
  return 0;
}

/* If startp and countp are not NULL they specify a
 * hyperslab of the variable to read.
 * startp gives the start index in each dimension
 * countp gives the number of entries in each dimension
 * If provided they must match with number_of_entries_per_dimension.
 *
 * pdata must point to an num_of_entries dimensional array of
 * unallocated data with size size_of_one_data_item.
 * The data will be allocated in this function.
 *
 * For example if number_of_dimensions = 2 and the type of data is double,
 * then pdata must be a double *** since it points to a double **.
 */
static int
t8_netcdf_read_data (const char *filename, const int ncid,
                     const char *varname,
                     const size_t number_of_dimensions,
                     const size_t * number_of_entries_per_dimension,
                     void **pdata, const size_t size_of_one_data_item,
                     const size_t * startp, const size_t * countp)
{
  int                 varid;
  int                 retval;
  int                 number_of_entries;
  size_t              idim;
  void               *data, **pdata_temp;

  /* Get the varid of the longitude data variable, based on its name. */
  t8_debugf ("Reading data info for '%s'\n", varname);
  retval = nc_inq_varid (ncid, varname, &varid);
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading data info", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  /* Ensure that the variable has the expected number of dimensions */
  {
    int                 ndims_in_file;
    retval = nc_inq_varndims (ncid, varid, &ndims_in_file);
    if (retval) {
      T8_NETCDF_ERROR (filename, "reading number of dimensions", retval);
    }
    if (ndims_in_file != (int) number_of_dimensions) {
      t8_global_errorf
        ("Error: '%s' variable does not have %zd dimension(s)\n",
         varname, number_of_dimensions);
      t8_netcdf_close_file (filename, ncid);
      return retval;
    }
    else {
      t8_debugf ("'%s' has exactly %zd dimension(s) as expected\n", varname,
                 number_of_dimensions);
    }
  }

#if 0
  /* TODO: - comment
   *       - move before file open
   *       - error handling
   */
  pdata_temp = pdata;
  number_of_entries = 1;
  for (idim = 0; idim < number_of_dimensions - 1; ++idim) {
    *pdata_temp =
      malloc (sizeof (void *) * number_of_entries_per_dimension[idim]);
    SC_CHECK_ABORT (*pdata_temp != NULL, "Could not allocate memory.");
    t8_debugf ("Allocated %zd void* on dim %zd\n",
               number_of_entries_per_dimension[idim], idim);
    pdata_temp = (void **) *pdata_temp;
    number_of_entries *= number_of_entries_per_dimension[idim];
  }
  *pdata_temp =
    malloc (size_of_one_data_item * number_of_entries_per_dimension[idim]);
  SC_CHECK_ABORT (*pdata_temp != NULL, "Could not allocate memory.");
  data = *pdata_temp;
  number_of_entries *= number_of_entries_per_dimension[idim];
  t8_debugf ("Allocated %zd data of size %zd on dim %zd\n",
             number_of_entries_per_dimension[idim], size_of_one_data_item,
             idim);

  if (data == NULL) {
    t8_global_errorf ("Could not allocate memory for %i data items\n",
                      number_of_entries);
    t8_netcdf_close_file (filename, ncid);
  }
  t8_debugf ("'%s' has %i entries\n", varname, number_of_entries);
#endif

  pdata_temp = pdata;
  number_of_entries = 1;
  for (idim = 0; idim < number_of_dimensions; ++idim) {
    number_of_entries *= number_of_entries_per_dimension[idim];
  }
  *pdata_temp = malloc (size_of_one_data_item * number_of_entries);
  data = *pdata_temp;
  if (data == NULL) {
    t8_global_errorf ("Could not allocate memory for %i data items\n",
                      number_of_entries);
    t8_netcdf_close_file (filename, ncid);
  }
  t8_debugf ("'%s' has %i entries\n", varname, number_of_entries);

  /* Read the data. */
  t8_debugf ("Reading '%s' data\n", varname);

  if (startp == NULL) {
    /* No hyperslab specified, read all data */
    T8_ASSERT (countp == NULL);
    t8_debugf ("reading all data\n");
    retval = nc_get_var (ncid, varid, data);
  }
  else {
    /* Read specified hyperslab of data */
    T8_ASSERT (startp != NULL && countp != NULL);
    t8_debugf ("reading hyperslab of data\n");
    retval = nc_get_vara (ncid, varid, startp, countp, data);
  }
  if (retval) {
    T8_NETCDF_ERROR (filename, "reading data", retval);
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }
  /* return success */
  return 0;
}

static inline       size_t
t8_netcdf_lat_long_to_point_index (const size_t lat_idx,
                                   const size_t long_idx,
                                   const size_t num_longitude_ixd)
{
  return num_longitude_ixd * lat_idx + long_idx;
}

#if 0                           /* Currently not needed */
static void
t8_netcdf_point_index_to_lat_long (const size_t point_idx,
                                   const size_t num_longitude_ixd,
                                   size_t * lat_idx_out,
                                   size_t * long_idx_out)
{
  *long_idx_out = point_idx % num_longitude_ixd;
  *lat_idx_out = point_idx / num_longitude_ixd;
}
#endif

static int
t8_netcdf_open_file (const char *filename, const double radius,
                     double **pcoordinates_euclidean,
                     size_t * pnum_latitude, size_t * pnum_longitude)
{
  int                 ncid, retval;
  int                 number_of_dims;
  int                 longitude_pos, latitude_pos;
  char                (*dimension_names)[BUFSIZ] = NULL;
  size_t             *dimension_lengths = NULL;
  float              *latitude_data, *longitude_data;
  int                *time;
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
    free (dimension_names);
    T8_FREE (dimension_lengths);
    return retval;
  }

  if (number_of_dims != 3) {
    t8_global_errorf ("Expected 3 dimensions. %i found.\n", number_of_dims);
    t8_netcdf_close_file (filename, ncid);
    free (dimension_names);
    T8_FREE (dimension_lengths);
    return 1;
  }

  /* Check the position of the longitude and latitude entries */
  longitude_pos = 0;
  latitude_pos = 1;
  T8_ASSERT (!strcmp (dimension_names[longitude_pos], "longitude"));
  T8_ASSERT (!strcmp (dimension_names[latitude_pos], "latitude"));
  T8_ASSERT (!strcmp (dimension_names[2], "time"));

  /* Read latitude data */
  retval =
    t8_netcdf_read_data (filename, ncid, "latitude", 1,
                         dimension_lengths + latitude_pos,
                         (void **) &latitude_data, sizeof (float), NULL,
                         NULL);
  if (retval) {
    /* An error occured and was printed,
     * the file is closed and we exit. */
    free (dimension_names);
    T8_FREE (dimension_lengths);
    return retval;
  }

  /* Read longitude data */
  retval =
    t8_netcdf_read_data (filename, ncid, "longitude", 1,
                         dimension_lengths + longitude_pos,
                         (void **) &longitude_data, sizeof (float), NULL,
                         NULL);
  if (retval) {
    /* An error occured and was printed,
     * the file is closed and we exit. */
    free (dimension_names);
    T8_FREE (dimension_lengths);
    return retval;
  }

  /* Read time data */
  retval =
    t8_netcdf_read_data (filename, ncid, "time", 1,
                         dimension_lengths + 2, (void **) &time,
                         sizeof (int), NULL, NULL);
  if (retval) {
    /* An error occured and was printed,
     * the file is closed and we exit. */
    free (dimension_names);
    T8_FREE (dimension_lengths);
    return retval;
  }

  /* Close the opened file */
  t8_netcdf_close_file (filename, ncid);

#ifdef T8_ENABLE_DEBUG
  /* Output the read data */
  {
    size_t              j;
    char                output[BUFSIZ] = "";
    char                number[20];
    for (j = 0; j < dimension_lengths[latitude_pos]; ++j) {
      snprintf (number, 20, " %.2f", latitude_data[j]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
    }
    t8_debugf ("Read data from '%s' with %zd entries:\n",
               dimension_names[latitude_pos],
               dimension_lengths[latitude_pos]);
    t8_debugf ("%s\n", output);

    strcpy (output, "");
    for (j = 0; j < dimension_lengths[longitude_pos]; ++j) {
      snprintf (number, 20, " %.2f", longitude_data[j]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
    }
    t8_debugf ("Read data from '%s' with %zd entries:\n",
               dimension_names[longitude_pos],
               dimension_lengths[longitude_pos]);
    t8_debugf ("%s\n", output);

    strcpy (output, "");
    for (j = 0; j < dimension_lengths[2]; ++j) {
      snprintf (number, 20, " %i", time[j]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
    }
    t8_debugf ("Read data from '%s' with %zd entries:\n",
               dimension_names[2], dimension_lengths[2]);
    t8_debugf ("%s\n", output);
  }
#endif

  /* Convert longitude and latitude to x,y,z */
  const size_t        num_long = dimension_lengths[longitude_pos];
  const size_t        num_lat = dimension_lengths[latitude_pos];
  *pnum_latitude = num_lat;
  *pnum_longitude = num_long;
  /* Compute the number of coordinates */
  const size_t        num_coordinates = num_long * num_lat;
  /* Allocate array to store all x,y,z coordinates */
  coordinates_euclidean = *pcoordinates_euclidean =
    T8_ALLOC (double, 3 * num_coordinates);

  /* Loop over all longitudes and all latitudes and compute the euclidean
   * coordinates for each point. */
  for (size_t ilat = 0; ilat < num_lat; ++ilat) {
    const double        latitude = latitude_data[ilat];
    for (size_t ilong = 0; ilong < num_long; ++ilong) {
      const double        longitude = longitude_data[ilong];
      double              xyz[3];
      /* Compute the current position in the euclidean array */
      const size_t        position =
        3 * t8_netcdf_lat_long_to_point_index (ilat, ilong, num_long);

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
  free (latitude_data);
  free (longitude_data);
  free (time);
  T8_FREE (dimension_lengths);
  free (dimension_names);

  /* Return success */
  return 0;
#undef NUM_DATA
}

/* Read msh-file and build a uniform forest on it.
 * Return 0 on success */
t8_forest_t
t8_reanalysis_build_forest (const char *mesh_filename, double radius,
                            int dimension, int level, sc_MPI_Comm comm)
{
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

  /* write the forest to vtk */
  {
    char                output_file_prefix[BUFSIZ];
    snprintf (output_file_prefix, BUFSIZ, "forest_uniform_l%i_%s",
              level, mesh_filename);
    t8_forest_write_vtk (forest, output_file_prefix);
  }

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
  const double       *coordinates;      /* The array of coordinates of all points. 
                                         * We assume that it is owned from outside and will not be destroyed if the search_data is destroyed. */
  sc_array_t         *matching_elements;        /* For each point an array of the element indices
                                                   that contain this point. (filled in the search query callback)
                                                 */
  sc_array_t         *points_per_element;       /* For each element an array of the point indices
                                                   of the points inside this elemnet (filled in the search query callback)
                                                 */
  size_t              num_points;       /* The number of points. coordinates has 3 * num_points entries
                                           and matching_elements has num_points entries. */
  t8_locidx_t         num_elements;     /* The number of elements. 
                                           points_per_element has num_elements entries. */
  int                 matched_elements; /* Count how many matching elements we find. */
  int                 matched_points;   /* Count for how many points we find at least one element. */
  double              point_finding_tolerance;  /* Numerical tolerance used for finding the points. 
                                                   The larger this value, the more points are found, 
                                                   even if they may not directly lie in the element, but are just close. */
  const double       *u10_data_per_element;     /* Interpolated u10 speed per element or NULL */
  const double       *v10_data_per_element;     /* Interpolated v10 speed per element or NULL */
  double              threshold;        /* Refinement threshold for refinement at u10v10 */
} t8_netcdf_search_data_t;

static void
t8_netcdf_destroy_search_data (t8_netcdf_search_data_t ** psearch_data)
{
  t8_netcdf_search_data_t *search_data;
  T8_ASSERT (psearch_data != NULL);
  search_data = *psearch_data;
  if (search_data != NULL) {
    size_t              ipoint;
    t8_locidx_t         ielement;

    /* Free the matching_elements array */
    for (ipoint = 0; ipoint < search_data->num_points; ++ipoint) {
      sc_array_reset (search_data->matching_elements + ipoint);
    }
    T8_FREE (search_data->matching_elements);

    /* Free the points_per_element array */
    for (ielement = 0; ielement < search_data->num_elements; ielement++) {
      sc_array_reset (search_data->points_per_element + ielement);
    }
    T8_FREE (search_data->points_per_element);
  }
  T8_FREE (search_data);
  *psearch_data = NULL;
}

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
  double             *tree_vertices;
  int                 is_definitely_outside;
  const double       *double_point = (const double *) point;

  /* Get a pointer to the vertex coordinates of the tree */
  tree_vertices = t8_forest_get_tree_vertices (forest, ltreeid);
  /* Do a quick estimate. This returns true if the point is definitely
   * outside of the element. */
  is_definitely_outside =
    t8_forest_element_point_outside_quick_estimate (forest, ltreeid, element,
                                                    tree_vertices,
                                                    double_point);

#if 1
  if (is_definitely_outside) {
    /* This point is not contained in this element, return 0 */
    return 0;
  }
#else
  is_definitely_outside = 0;
#endif

  t8_netcdf_search_data_t *user_data =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest);
  /* The point may be inside the element. Do a proper check. */
  if (t8_forest_element_point_inside
      (forest, ltreeid, element, double_point,
       user_data->point_finding_tolerance)) {
    T8_ASSERT (!is_definitely_outside); /* Should never happen */
    /* This point is contained in this element */
    if (is_leaf) {
      /* This element is a leaf element, we add its index to the list of
       * elements that contain this point */
      /* In debugging mode we count how many points and elements we match */
      user_data->matched_elements++;
      if (user_data->matching_elements[point_index].elem_count == 0) {
        /* This point was not found inside an element yet, we add to the counter of matched points */
        user_data->matched_points++;
      }

      /* Compute the forest local index of the element */
      t8_locidx_t         element_index =
        t8_forest_get_tree_element_offset (forest, ltreeid) + tree_leaf_index;

      t8_debugf ("Pushing index %i to point %zd\n", element_index,
                 point_index);
      /* Add this index to the array of found elements */
      *(t8_locidx_t *) sc_array_push (user_data->matching_elements +
                                      point_index) = element_index;
      /* Add the point's index to the array of points matching this element */
      *(size_t *) sc_array_push (user_data->points_per_element +
                                 element_index) = point_index;
    }
    /* Since the point is contained in the element, we return 1 */
    return 1;
  }
  else {
    /* This point is not contained in this element, return 0 */
    return 0;
  }
}

static void
t8_netcdf_init_search_data (t8_forest_t forest, const double *points,
                            const size_t num_points, const double tolerance,
                            t8_netcdf_search_data_t * search_data)
{
  sc_array_t         *matching_elements;
  sc_array_t         *points_per_elements;
  t8_locidx_t         num_elements, ielement;
  size_t              ipoint;

  /* Allocate as many arrays as we have points to store the
   * matching elements for each point */
  matching_elements = T8_ALLOC (sc_array_t, num_points);
  /* Initialize these arrays */
  for (ipoint = 0; ipoint < num_points; ++ipoint) {
    sc_array_init (matching_elements + ipoint, sizeof (t8_locidx_t));
  }
  /* Allocate the points per element arrays */
  num_elements = t8_forest_get_local_num_elements (forest);
  points_per_elements = T8_ALLOC (sc_array_t, num_elements);
  /* Initialize these arrays */
  for (ielement = 0; ielement < num_elements; ++ielement) {
    sc_array_init (points_per_elements + ielement, sizeof (size_t));
  }

  /* Init user data for the search routine */
  search_data->num_points = num_points;
  search_data->num_elements = num_elements;
  search_data->coordinates = points;
  search_data->matching_elements = matching_elements;
  search_data->points_per_element = points_per_elements;
  search_data->point_finding_tolerance = tolerance;
  /* Set start/default values */
  search_data->matched_elements = 0;
  search_data->matched_points = 0;
  search_data->u10_data_per_element = NULL;
  search_data->v10_data_per_element = NULL;
  search_data->threshold = 0;
}

/* Given a forest and an array of points, identify the elements that contain
 * the points.
 * The points are given as one coordinate array in the format (x_0 y_0 z_0 x_1 y_1 z_1 ... )
 */
t8_netcdf_search_data_t *
t8_netcdf_find_mesh_elements (t8_forest_t forest, double *points,
                              const size_t num_points, const double tolerance)
{
  t8_netcdf_search_data_t *coords_and_matching_elements;
  sc_array_t          queries;
  double              search_runtime;

  /* Allocate search data */
  coords_and_matching_elements = T8_ALLOC (t8_netcdf_search_data_t, 1);

  /* Initialize search data */
  t8_netcdf_init_search_data (forest, points, num_points, tolerance,
                              coords_and_matching_elements);

  /* Set this data as the forests user pointer */
  t8_forest_set_user_data (forest, coords_and_matching_elements);

  /* Initialize the array of points to be passed to the search function.
   * Each entry is one point, thus 3 doubles */
  sc_array_init_data (&queries, points, 3 * sizeof (double), num_points);
  t8_global_productionf ("Starting search with %zd points\n", num_points);
  search_runtime = -sc_MPI_Wtime ();
  t8_forest_search (forest, t8_netcdf_find_mesh_elements_query,
                    t8_netcdf_find_mesh_elements_query, &queries);
  search_runtime += sc_MPI_Wtime ();
  sc_array_reset (&queries);

  t8_global_productionf
    ("Finished search. Found %i points and matched %i elements\n",
     coords_and_matching_elements->matched_points,
     coords_and_matching_elements->matched_elements);
  t8_global_productionf ("Search runtime: %.3fs\n", search_runtime);

  return coords_and_matching_elements;
}

typedef struct
{
  double             *u10_data_per_element;
  double             *v10_data_per_element;
} t8_netcdf_u10v10data_t;

/* Refine callback that refines an element if its combined u10 and v10
 * speed is larger (in absolute value) than a threshold. */
static              t8_locidx_t
t8_netcdf_refine_if_u10v10_large (t8_forest_t forest,
                                  t8_forest_t forest_from,
                                  t8_locidx_t which_tree,
                                  t8_locidx_t lelement_id,
                                  t8_eclass_scheme_c * ts,
                                  int num_elements, t8_element_t * elements[])
{
  const t8_netcdf_search_data_t *u10v10data =
    (const t8_netcdf_search_data_t *) t8_forest_get_user_data (forest);
  t8_locidx_t         element_index;
  const double        threshold = u10v10data->threshold;        /* Refine if absolute value is larger */
  double              absolute_val;

  /* Compute the index of the current element */
  element_index =
    t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

  absolute_val = fabs (u10v10data->u10_data_per_element[element_index])
    + fabs (u10v10data->v10_data_per_element[element_index]);

  if (absolute_val > threshold) {
    /* Refine this element */
    return 1;
  }
  return 0;
}

static              t8_locidx_t
t8_netcdf_count_double_matched_elements (const t8_netcdf_search_data_t *
                                         search_result)
{
  t8_locidx_t         ielement, num_double_matched = 0;

  for (ielement = 0; ielement < search_result->num_elements; ++ielement) {
    if (search_result->points_per_element[ielement].elem_count > 1) {
      num_double_matched++;
    }
  }
  return num_double_matched;
}

/* Refine callback that refines an element if multiple points were found
 * for it. */
static              t8_locidx_t
t8_netcdf_refine_if_multiple_match (t8_forest_t forest,
                                    t8_forest_t forest_from,
                                    t8_locidx_t which_tree,
                                    t8_locidx_t lelement_id,
                                    t8_eclass_scheme_c * ts,
                                    int num_elements,
                                    t8_element_t * elements[])
{
  t8_netcdf_search_data_t *search_result;
  t8_locidx_t         element_index;

  /* Get the result of the last search */
  search_result =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest_from);
  /* Compute the index of the current element */
  element_index =
    t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

  if (search_result->points_per_element[element_index].elem_count > 1) {
    /* This element matched more than one point, we refine it */
    return 1;
  }
  return 0;
}

/* Refine the forest in the elements that have multiple points found */
static              t8_forest_t
t8_netcdf_refine (t8_forest_t forest)
{
  t8_forest_t         forest_refine =
    t8_forest_new_adapt (forest, t8_netcdf_refine_if_multiple_match,
                         0, 0, NULL);

  return forest_refine;
}

static void
t8_netcdf_interpolate_to_refined_forest (t8_forest_t forest_old,
                                         t8_forest_t forest_new,
                                         t8_locidx_t which_tree,
                                         t8_eclass_scheme_c * ts,
                                         int num_outgoing,
                                         t8_locidx_t first_outgoing,
                                         int num_incoming,
                                         t8_locidx_t first_incoming)
{
  t8_netcdf_search_data_t *search_data_old;
  t8_netcdf_search_data_t *search_data_new;
  t8_locidx_t         element_index_old, element_index_new;
  sc_array_t         *points_per_element_old;
  sc_array_t         *points_per_element_new;
  t8_locidx_t         ichild;

  search_data_old =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest_old);
  search_data_new =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest_new);

  element_index_old =
    t8_forest_get_tree_element_offset (forest_old,
                                       which_tree) + first_outgoing;
  element_index_new =
    t8_forest_get_tree_element_offset (forest_new,
                                       which_tree) + first_incoming;

  /* Get a pointer to the points per element array of the old and new element */
  points_per_element_old =
    &search_data_old->points_per_element[element_index_old];
  points_per_element_new =
    &search_data_new->points_per_element[element_index_new];

  size_t              ipoint, point_index;
  const double       *point;
  t8_element_t       *child, *element;

  /* Loop over all points */
  for (ipoint = 0; ipoint < points_per_element_old->elem_count; ++ipoint) {
    /* Get the position of the point and get its coordinates */
    point_index = *(size_t *) sc_array_index (points_per_element_old, ipoint);
    point = &search_data_old->coordinates[3 * point_index];
    t8_debugf ("Searching point %zd in children of elem %i\n", point_index,
               element_index_old);

    if (num_incoming == num_outgoing) {
      T8_ASSERT (num_incoming == 1);
      /* This element was not refined, we search again for
       * the points that we initially found, but with the new tolerande. */
      element =
        t8_forest_get_element_in_tree (forest_new, which_tree,
                                       first_incoming);
      if (t8_forest_element_point_inside
          (forest_new, which_tree, element, point,
           search_data_new->point_finding_tolerance)) {
        /* This point is in this element, we add it to the new list */
        *(size_t *) sc_array_push (points_per_element_new) = point_index;
      }
    }
    else if (num_incoming > 1) {
      /* This element is refined into num_incoming children */
      /* For each point and each child we check whether the point is within the child. */

      /* Loop over all children */
      for (ichild = 0; ichild < num_incoming; ++ichild) {
        t8_debugf ("child %i\n", ichild);
        child =
          t8_forest_get_element_in_tree (forest_new, which_tree,
                                         first_incoming + ichild);
        if (t8_forest_element_point_inside
            (forest_new, which_tree, child, point,
             search_data_new->point_finding_tolerance)) {
          t8_debugf ("found on child %i\n", ichild);
          /* This point is in this element, we add it to the new list */
          *(size_t *) sc_array_push (points_per_element_new + ichild) =
            point_index;
        }
      }
    }
#if T8_ENABLE_DEBUG
    else {
      /* The element is coarsened. This should not happen. */
      T8_ASSERT (0 && "Element is coarsened. This should not happen");
    }
#endif
  }
  for (ichild = 0; ichild < num_incoming; ++ichild) {
    if ((points_per_element_new + ichild)->elem_count == 0
        && points_per_element_old->elem_count > 0) {
      /* There were points on this element before, but we did not find any now.
       * In this case, we copy all points of the old element. */
      sc_array_copy (points_per_element_new + ichild, points_per_element_old);
    }
  }
}

/* Given the original forest and a refined forest, we allocate
 * new search data for the refined forest and interpolate the 
 * data from the non refined forest.
 * In particular for the elements that matched multiple points in the
 * first forest we refined them and try to map the resulting points to
 * the newly refined elements. */
static void
t8_netcdf_interpolate_data (t8_forest_t forest_refined,
                            t8_forest_t forest_not_refined)
{
  t8_netcdf_search_data_t *search_data_old;
  t8_netcdf_search_data_t *search_data_new;

  search_data_old =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest_not_refined);
  /* Allocate new search data */
  search_data_new = T8_ALLOC (t8_netcdf_search_data_t, 1);
  /* Initialize new search data.
   * We use the coordinates array of the old search data here. 
   * We also use a smaller numerical tolerance. */
  t8_netcdf_init_search_data (forest_refined, search_data_old->coordinates,
                              search_data_old->num_points,
                              search_data_old->point_finding_tolerance * .7,
                              search_data_new);

  /* Set this new search data as the user pointer of the refined forest */
  t8_forest_set_user_data (forest_refined, search_data_new);

  t8_forest_iterate_replace (forest_refined, forest_not_refined,
                             t8_netcdf_interpolate_to_refined_forest);

}

static int
t8_netcdf_read_data_to_forest (const char *filename, t8_forest_t forest,
                               const size_t num_latitude,
                               const size_t num_longitude,
                               const size_t num_timesteps,
                               t8_netcdf_search_data_t *
                               search_results, sc_MPI_Comm comm)
{
  double             *u10_data_per_element;
  double             *v10_data_per_element;
  short              *u10_data_from_file;
  short              *v10_data_from_file;
  int                 retval;
  int                 ncid, u10varid, v10varid;
  const size_t        startp[3] = { 0, 0, 0 };
  const size_t        countp[3] =
    { num_timesteps, num_latitude, num_longitude };
  const size_t        num_data_per_dimension[3] =
    { num_timesteps, num_latitude, num_longitude };
  const size_t        num_points = num_longitude * num_latitude;
  t8_locidx_t         num_elements, matched_elements, global_matched_elements;
  double              u10_scaling, v10_scaling;
  double              u10_offset, v10_offset;
  size_t              time_step;
  t8_locidx_t         ielement;
  size_t              ipoint;

  T8_ASSERT (num_points <=
             num_data_per_dimension[0] * num_data_per_dimension[1]
             * num_data_per_dimension[2]);

  /* Open the file */
  t8_debugf ("Opening file %s\n", filename);
  retval = nc_open (filename, NC_NOWRITE, &ncid);
  if (retval) {
    /* Could not open the file */
    T8_NETCDF_ERROR (filename, "opening file", retval);
    return retval;
  }

  /* Get the varid of u10 */
  retval = nc_inq_varid (ncid, "u10", &u10varid);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading variable id of u10", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  retval =
    t8_netcdf_read_data (filename, ncid, "u10", 3, num_data_per_dimension,
                         (void **) &u10_data_from_file,
                         sizeof (*u10_data_from_file), startp, countp);
  if (retval) {
    /* Could not read the data.
     * The file was closed in the error handling of t8_netcdf_read_data. */
    T8_NETCDF_ERROR (filename, "reading u10 data", retval);
    return retval;
  }

  /* Read the scaling and offset attributes */
  retval = nc_get_att (ncid, u10varid, "scale_factor", &u10_scaling);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading u10:scale_factor attribute", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }
  /* offset */
  retval = nc_get_att (ncid, u10varid, "add_offset", &u10_offset);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading u10:add_offset attribute", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  /* Get the varid of v10 */
  retval = nc_inq_varid (ncid, "v10", &v10varid);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading variable id of v10", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  retval =
    t8_netcdf_read_data (filename, ncid, "v10", 3, num_data_per_dimension,
                         (void **) &v10_data_from_file,
                         sizeof (*v10_data_from_file), startp, countp);
  if (retval) {
    /* Could not read the data.
     * The file was closed in the error handling of t8_netcdf_read_data. */
    T8_NETCDF_ERROR (filename, "reading u10 data", retval);
    return retval;
  }

  /* Read the scaling and offset attributes */
  retval = nc_get_att (ncid, v10varid, "scale_factor", &v10_scaling);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading v10:scale_factor attribute", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }
  /* offset */
  retval = nc_get_att (ncid, v10varid, "add_offset", &v10_offset);
  if (retval) {
    /* Error */
    T8_NETCDF_ERROR (filename, "reading v10:add_offset attribute", retval);
    /* close the file */
    t8_netcdf_close_file (filename, ncid);
    return retval;
  }

  /* We are done reading from the file and thus close it */
  t8_netcdf_close_file (filename, ncid);

  t8_debugf ("Read u10 data with scaling %g and offset %g\n", u10_scaling,
             u10_offset);
  t8_debugf ("Read v10 data with scaling %g and offset %g\n", v10_scaling,
             v10_offset);
#if 1
//#ifdef T8_ENABLE_DEBUG
  size_t              num_print_data = 100;     // only print the first N values
  /* Output the read data */
  {
    size_t              j;
    char                output[BUFSIZ] = "";
    char                number[20];
    t8_debugf ("Read data from 'u10' with %zd entries:\n", num_points);
    for (j = 0; j < num_print_data; ++j) {
      snprintf (number, 20, " %i", u10_data_from_file[j]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
      else {
        t8_debugf ("%s\n", output);
        /* Overwrite output */
        strcpy (output, number);
      }
    }
    t8_debugf ("%s\n", output);
    output[0] = '\0';
    t8_debugf ("Read data from 'v10' with %zd entries:\n", num_points);
    for (j = 0; j < num_print_data; ++j) {
      snprintf (number, 20, " %i", v10_data_from_file[j]);
      if (strlen (output) < BUFSIZ - 21) {
        strcat (output, number);
      }
      else {
        t8_debugf ("%s\n", output);
        /* Overwrite output */
        strcpy (output, number);
      }
    }
    t8_debugf ("%s\n", output);
  }
//#endif
#endif

  /* We now fill the per element array with the data */
  num_elements = t8_forest_get_local_num_elements (forest);
  u10_data_per_element = T8_ALLOC_ZERO (double, num_elements);
  v10_data_per_element = T8_ALLOC_ZERO (double, num_elements);

  t8_debugf ("Projecting data to %i elements\n", num_elements);

  time_step = 0;
  matched_elements = 0;
  T8_ASSERT (num_elements == search_results->num_elements);
  for (ielement = 0; ielement < num_elements; ++ielement) {
    if (search_results->points_per_element[ielement].elem_count > 0) {
      /* We found a point for this element */
      /* We average the values of all these points and add them to the data. */
      size_t              num_points_in_element =
        search_results->points_per_element[ielement].elem_count;
      for (ipoint = 0; ipoint < num_points_in_element; ++ipoint) {
        size_t              point_idx = *(size_t *)
          sc_array_index (&search_results->points_per_element[ielement],
                          ipoint);

        /* Add all these values up */
        u10_data_per_element[ielement] +=
          u10_data_from_file[point_idx] * u10_scaling + u10_offset;
        v10_data_per_element[ielement] +=
          v10_data_from_file[point_idx] * v10_scaling + v10_offset;
      }
      /* Divdie by the number of values */
      u10_data_per_element[ielement] /= num_points_in_element;
      v10_data_per_element[ielement] /= num_points_in_element;
      T8_ASSERT (time_step == 0);       /* Currently time_step = 0 is fixed */
      matched_elements++;       /* Count this element as matched */

      t8_debugf ("Writing %g to element %i from %zd points of %zd\n",
                 u10_data_per_element[ielement],
                 ielement, num_points_in_element, num_points);
    }
    else {
      /* There is no point on this element. We store a large negative value in v10,
         such that we can identify the element in the vtk. */
      v10_data_per_element[ielement] = -10000;
    }
  }

  sc_MPI_Reduce (&matched_elements, &global_matched_elements, 1, sc_MPI_INT,
                 sc_MPI_SUM, 0, comm);
  t8_global_productionf ("Interpolated to %i elements.\n",
                         global_matched_elements);

  /* Store the interpolated u10 and v10 values in the search result */
  search_results->u10_data_per_element = u10_data_per_element;
  search_results->v10_data_per_element = v10_data_per_element;

  t8_global_productionf ("Stored at %p\n", (void *) search_results);
  t8_global_productionf ("At 100: %g  %g\n", u10_data_per_element[100],
                         v10_data_per_element[100]);

#if 0
  /* Refine */
  t8_forest_t         forest_refine;
  {
    t8_netcdf_u10v10data_t refine_data =
      { u10_data_per_element, v10_data_per_element };

    t8_forest_ref (forest);
    forest_refine =
      t8_forest_new_adapt (forest, t8_netcdf_refine_if_u10v10_large, 0, 0,
                           &refine_data);

    t8_netcdf_interpolate_data (forest_refine, forest);
    //t8_forest_write_vtk (forest_refine, "test_u10v10_refined");
    //t8_forest_unref (&forest_refine);
    forest = forest_refine;
  }
#endif

  /* clean-up */
  free (u10_data_from_file);
  free (v10_data_from_file);

  return 0;
}

    /* VTK output
     * We assume that a search_data object is set as the forest's user data */
static void
t8_netcdf_write_vtk (t8_forest_t forest, const char *fileprefix)
{
  t8_netcdf_search_data_t *search_result =
    (t8_netcdf_search_data_t *) t8_forest_get_user_data (forest);
  t8_vtk_data_field_t vtk_data[2];
  vtk_data[0].data = (double *) search_result->u10_data_per_element;
  vtk_data[1].data = (double *) search_result->v10_data_per_element;
  T8_ASSERT (vtk_data[0].data != NULL && vtk_data[1].data != NULL);
  strcpy (vtk_data[0].description, "u10");
  strcpy (vtk_data[1].description, "v10");
  vtk_data[0].type = T8_VTK_SCALAR;
  vtk_data[1].type = T8_VTK_SCALAR;
  t8_forest_vtk_write_file (forest, fileprefix, 1, 1, 1, 1, 0, 2, vtk_data);
}

static void
t8_netcdf_adapt_spheres (const char *mesh_filename, const int sphere_radius,
                         const int sphere_dim, const int level,
                         const int refine_levels, sc_MPI_Comm comm)
{
  /* Build uniform forest from mesh file */
  t8_forest_t         forest =
    t8_reanalysis_build_forest (mesh_filename, sphere_radius, sphere_dim,
                                level,
                                comm);
  t8_forest_t         forest_adapt = NULL;
  int                 ref_level;

  /* Define Spheres to refine along */
  int                 num_spheres = 8;  /* Number of spheres */
  /* Initialize spheres array */
  t8_example_sphere_t *sphere;
  sc_array_t         *spheres =
    sc_array_new_count (sizeof (t8_example_sphere_t), num_spheres);
  /* Set sphere 1 (island) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 0);
  sphere->midpoint[0] = 0.225;
  sphere->midpoint[1] = -0.07;
  sphere->midpoint[2] = 0.45;
  sphere->radius = 0.06;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 2 (close to island) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 1);
  sphere->midpoint[0] = 0.245;
  sphere->midpoint[1] = -0.07;
  sphere->midpoint[2] = 0.435;
  sphere->radius = 0.06;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 3 (away to the left 1) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 2);
  sphere->midpoint[0] = 0.255;
  sphere->midpoint[1] = -0.1;
  sphere->midpoint[2] = 0.42;
  sphere->radius = 0.04;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 4 (away to the left 2) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 3);
  sphere->midpoint[0] = 0.26;
  sphere->midpoint[1] = -0.115;
  sphere->midpoint[2] = 0.415;
  sphere->radius = 0.05;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 5 (away to the left 3) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 4);
  sphere->midpoint[0] = 0.27;
  sphere->midpoint[1] = -0.15;
  sphere->midpoint[2] = 0.4;
  sphere->radius = 0.06;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 6 (coarse ocean) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 5);
  sphere->midpoint[0] = 0.33;
  sphere->midpoint[1] = -0.08;
  sphere->midpoint[2] = 0.39;
  sphere->radius = 0.095;
  sphere->max_level = level + refine_levels - 1;
  sphere->refine_outside = 0;

  /* Set sphere 7 (France 1) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 6);
  sphere->midpoint[0] = 0.33;
  sphere->midpoint[1] = 0.01;
  sphere->midpoint[2] = 0.38;
  sphere->radius = 0.05;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  /* Set sphere 8 (France 2) */
  sphere = (t8_example_sphere_t *) sc_array_index (spheres, 7);
  sphere->midpoint[0] = 0.34;
  sphere->midpoint[1] = 0.03;
  sphere->midpoint[2] = 0.37;
  sphere->radius = 0.07;
  sphere->max_level = level + refine_levels;
  sphere->refine_outside = 0;

  for (ref_level = level; ref_level < refine_levels + level; ++ref_level) {
    forest_adapt =
      t8_forest_new_adapt (forest, t8_common_adapt_union_of_spheres, 0, 0,
                           spheres);

    /* Scale the spheres */
    size_t              isphere;
    for (isphere = 0; isphere < spheres->elem_count; ++isphere) {
      t8_example_sphere_t *sphere =
        (t8_example_sphere_t *) sc_array_index (spheres, isphere);
      sphere->radius *= ref_level == level ? 0.9 : 0.5;
    }
    forest = forest_adapt;
  }

  t8_forest_init (&forest);
  t8_forest_set_balance (forest, forest_adapt, 0);
  t8_forest_commit (forest);

  /* write vtk */
  t8_forest_write_vtk (forest, "adapt_spheres");
  /* clean up */
  sc_array_destroy (spheres);
  t8_forest_unref (&forest);
}

static void
t8_netcdf_read_data_and_adapt (const char *netcdf_filename,
                               const char *mesh_filename,
                               const double sphere_radius,
                               const int sphere_dim, const int level,
                               sc_MPI_Comm comm)
{
  int                 retval;
  t8_forest_t         forest =
    t8_reanalysis_build_forest (mesh_filename, sphere_radius, sphere_dim,
                                level,
                                comm);
  t8_forest_t         forest_refined;

  if (forest != NULL) {
    double             *coordinates_euclidean;
    size_t              num_coordinates;
    size_t              num_latitude, num_longitude;
    t8_netcdf_search_data_t *search_result;
    const double        tolerance = 1e-1;

    retval = t8_netcdf_open_file (netcdf_filename, sphere_radius,
                                  &coordinates_euclidean, &num_latitude,
                                  &num_longitude);
    num_coordinates = num_latitude * num_longitude;
    if (!retval) {
#if 1
#ifdef T8_ENABLE_DEBUG
      size_t              max_num_coordinates = num_latitude;
      size_t              new_num_coordinates =
        SC_MIN (num_coordinates, max_num_coordinates);
      t8_debugf
        ("Debugging mode detected. Search only for %zd of the %zd points"
         " in order to reduce the runtime (debugging mode is slow).\n",
         new_num_coordinates, num_coordinates);
      num_coordinates = new_num_coordinates;
      if (new_num_coordinates == max_num_coordinates) {
        num_longitude = 1;
      }
#endif
#endif
      search_result =
        t8_netcdf_find_mesh_elements (forest, coordinates_euclidean,
                                      num_coordinates, tolerance);

      /* Print double matched elements */
      t8_global_productionf ("Matched %i elements at least twice\n",
                             t8_netcdf_count_double_matched_elements
                             (search_result));
      /* We now refine the forest */
      t8_forest_ref (forest);   /* Ref the original forest since we need to use it after refinement */
      /* Refine the forest in the elements that match multiple points */
      forest_refined = t8_netcdf_refine (forest);
      t8_netcdf_interpolate_data (forest_refined, forest);
      t8_forest_write_vtk (forest, "test_prerefine");
      /* Print double matched elements */
      search_result = (t8_netcdf_search_data_t *)
        t8_forest_get_user_data (forest_refined);
      t8_global_productionf ("Matched %i elements at least twice\n",
                             t8_netcdf_count_double_matched_elements
                             (search_result));
#if 0
      /* Do a second run, refine/interpolate */
      t8_forest_unref (&forest);
      forest = forest_refined;
      t8_forest_ref (forest);
      forest_refined = t8_netcdf_refine (forest);
      t8_netcdf_interpolate_data (forest_refined, forest);
      t8_forest_write_vtk (forest, "test_prerefine2");
      /* Print double matched elements */
      search_result = (t8_netcdf_search_data_t *)
        t8_forest_get_user_data (forest_refined);
      t8_global_productionf ("Matched %i elements at least twice\n",
                             t8_netcdf_count_double_matched_elements
                             (search_result));
#endif
      /* Do a third run */
      t8_forest_unref (&forest);
      forest = forest_refined;
      t8_forest_ref (forest);
      forest_refined = t8_netcdf_refine (forest);
      t8_netcdf_interpolate_data (forest_refined, forest);
      t8_forest_write_vtk (forest, "test_prerefine3");

      /* The old search_result was destroyed in the interpolate routine.
       * Instead the updated search result is stored at the refined forest. */
      search_result = (t8_netcdf_search_data_t *)
        t8_forest_get_user_data (forest_refined);

      t8_netcdf_read_data_to_forest (netcdf_filename, forest_refined,
                                     num_latitude, num_longitude, 1,
                                     search_result, comm);

      t8_global_productionf ("Matched %i elements at least twice\n",
                             t8_netcdf_count_double_matched_elements
                             (search_result));

      t8_netcdf_write_vtk (forest_refined, "test");

      /* Refine the forest where u10 and v10 are large */
      t8_global_productionf ("Start refining at u10v10\n");
      /* Set threshold for refinement */
      search_result->threshold = 6;
      t8_forest_ref (forest_refined);
      t8_forest_t         forest_adapt;
      forest_adapt =
        t8_forest_new_adapt (forest_refined,
                             t8_netcdf_refine_if_u10v10_large,
                             0, 0, search_result);

      t8_netcdf_interpolate_data (forest_adapt, forest_refined);
      /* The old search_result was destroyed in the interpolate routine.
       * Instead the updated search result is stored at the refined forest. */
      search_result = (t8_netcdf_search_data_t *)
        t8_forest_get_user_data (forest_adapt);
      /* Clean up the old element data */
      T8_FREE ((void *) search_result->u10_data_per_element);
      T8_FREE ((void *) search_result->v10_data_per_element);

      t8_netcdf_read_data_to_forest (netcdf_filename, forest_adapt,
                                     num_latitude, num_longitude, 1,
                                     search_result, comm);
      /* Do another run */
      t8_forest_ref (forest_refined);
      t8_forest_ref (forest_adapt);
      t8_forest_t         forest_adapt2;
      /* Set threshold for refinement */
      search_result->threshold = 15;
      forest_adapt2 =
        t8_forest_new_adapt (forest_adapt,
                             t8_netcdf_refine_if_u10v10_large,
                             0, 0, search_result);

      t8_netcdf_interpolate_data (forest_adapt2, forest_adapt);
      /* The old search_result was destroyed in the interpolate routine.
       * Instead the updated search result is stored at the refined forest. */
      search_result = (t8_netcdf_search_data_t *)
        t8_forest_get_user_data (forest_adapt2);
      /* Clean up the old element data */
      T8_FREE ((void *) search_result->u10_data_per_element);
      T8_FREE ((void *) search_result->v10_data_per_element);

      t8_netcdf_read_data_to_forest (netcdf_filename, forest_adapt2,
                                     num_latitude, num_longitude, 1,
                                     search_result, comm);
      /* Repartition */
      t8_forest_ref (forest_adapt2);
      t8_forest_t         forest_partition;
      sc_array_t          u10data_old_view;
      sc_array_t          v10data_old_view;
      sc_array_t          u10data_new_view;
      sc_array_t          v10data_new_view;
      t8_locidx_t         num_elements_partitioned;
      double             *u10data_new, *v10data_new;
      sc_array_init_data (&u10data_old_view,
                          (void *) search_result->u10_data_per_element,
                          sizeof (*search_result->u10_data_per_element),
                          search_result->num_elements);
      sc_array_init_data (&v10data_old_view,
                          (void *) search_result->u10_data_per_element,
                          sizeof (*search_result->u10_data_per_element),
                          search_result->num_elements);
      t8_forest_init (&forest_partition);
      t8_forest_set_partition (forest_partition, forest_adapt2, 0);
      t8_forest_commit (forest_partition);
      num_elements_partitioned =
        t8_forest_get_local_num_elements (forest_partition);
      u10data_new = T8_ALLOC (double, num_elements_partitioned);
      v10data_new = T8_ALLOC (double, num_elements_partitioned);
      sc_array_init_data (&u10data_new_view, u10data_new, sizeof (double),
                          num_elements_partitioned);
      sc_array_init_data (&v10data_new_view, v10data_new, sizeof (double),
                          num_elements_partitioned);
      t8_forest_partition_data (forest_adapt2, forest_partition,
                                &u10data_old_view, &u10data_new_view);
      t8_forest_partition_data (forest_adapt2, forest_partition,
                                &v10data_old_view, &v10data_new_view);
      search_result->u10_data_per_element = u10data_new;
      search_result->v10_data_per_element = v10data_new;
      /* TODO: Careful, search_result does now only partly match the new forest.
         the num_elements and matched_elements etc. arrays are still with the old forest. */
      t8_forest_set_user_data (forest_partition, search_result);
      t8_netcdf_write_vtk (forest_partition,
                           "test_refinedu10v10_partitioned");

      /* TODO: Store indices instead of values? interpolate must interpolated values to? */

      /* Clean-up */
      T8_FREE (coordinates_euclidean);
      T8_FREE ((void *) search_result->u10_data_per_element);
      T8_FREE ((void *) search_result->v10_data_per_element);
      t8_netcdf_destroy_search_data (&search_result);
      t8_forest_unref (&forest_refined);
      t8_forest_unref (&forest_adapt);
      t8_forest_unref (&forest_adapt2);
      t8_forest_unref (&forest_partition);
    }
    t8_forest_unref (&forest);
  }
}
#endif

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  const char         *netcdf_filename = NULL;
  const char         *mesh_filename = NULL;
  int                 parsed, helpme;
  int                 sphere_dim, level;
  int                 adapt_spheres;
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
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The uniform refinement level of the mesh. Default = 0");
  sc_options_add_int (opt, 's', "s", &adapt_spheres, 0,
                      "Adapt the given number of levels according to given sphere coordinates (experimental). Deactivates data interpolation. Default = 0");

  /* Parse the command line arguments from the input */
  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);

  if (parsed >= 0 && helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
#if T8_WITH_NETCDF
  else if (parsed >= 0 && netcdf_filename != NULL && mesh_filename != NULL
           && 0 <= level && 0 <= adapt_spheres) {
    if (!adapt_spheres) {
      t8_netcdf_read_data_and_adapt (netcdf_filename, mesh_filename,
                                     sphere_radius, sphere_dim, level, comm);
    }
    else {
      t8_netcdf_adapt_spheres (mesh_filename, sphere_radius, sphere_dim,
                               level, adapt_spheres, comm);
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
