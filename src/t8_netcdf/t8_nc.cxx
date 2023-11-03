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

#include "t8_netcdf/t8_nc.h"
#include <t8_netcdf/t8_nc_data.hxx>
#include <t8_netcdf/t8_nc_mesh.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>

/* An enumerator for array indices of at most four dimensional geo-spatial (and temporal) data */
enum t8_coord_ids { T8_COORD_ID_UNDEFINED = -1, T8_LON = 0, T8_LAT = 1, T8_LEV = 2, T8_TIME = 3, T8_NUM_COORD_IDS = 4 };

/* A macro for an internal error corresponding to netCDF functionality */
#define T8_NC_ERR_CODE -1

/* Define a default value for not considered information */
#define T8_NC_NOT_CONSIDERED -1

struct t8_nc_data
{
 private:
  const int ncid;

 public:
  t8_nc_data (const int _ncid): ncid { _ncid } {};
  t8_nc_data (const int _ncid, const bool parallel_access): ncid { _ncid }, use_distributed_data { parallel_access } {};
  t8_nc_data (const int _ncid, const bool parallel_access, const sc_MPI_Comm _comm)
    : ncid { _ncid }, use_distributed_data { parallel_access }, comm { _comm } {};
  ~t8_nc_data ()
  {
    for (auto iter { coordinates.begin () }; iter != coordinates.end (); ++iter) {
      /* Deallocate the coordinate variable */
      if (*iter != nullptr) {
        t8_nc_destroy_geo_variable (*iter);
      }
    }
    for (auto iter { vars.begin () }; iter != vars.end (); ++iter) {
      /* Deallocate the variable */
      t8_nc_destroy_geo_variable (*iter);
    }
    if (mesh != nullptr) {
      /* Deallocate the constructed mesh */
      t8_nc_mesh_destroy (mesh);
    }
  };

  /* Save the geo-spatial coordinate IDs */
  std::array<int, t8_coord_ids::T8_NUM_COORD_IDS> coord_dim_ids {
    T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED
  };  //!< An array storing the dimension IDs of the geo-spatial dimensions
  std::array<int, t8_coord_ids::T8_NUM_COORD_IDS> coord_var_ids {
    T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED, T8_NC_NOT_CONSIDERED
  };  //!< An array storing the variable IDs of the coordinate variables (compliant to the geo-spatial dimensions)
  std::array<size_t, t8_coord_ids::T8_NUM_COORD_IDS> coord_lengths {
    1, 1, 1, 1
  };  //!< An array storing the length of each geo-spatial dimension
  std::array<t8_geo_var_t, t8_coord_ids::T8_NUM_COORD_IDS> coordinates {
    nullptr, nullptr, nullptr, nullptr
  };  //!< An array holding the actual coordinate data of the coordinate variables

  /* A map for storing hints regarding the dimension (i.e. their names) */
  std::map<std::basic_string<char>, t8_coord_ids>
    dimension_hints;  //!< A vector storing user-specified hints about the geo-spatial dimension (Advised if the dimensions do not have standard names like 'longitude' etc.)

  /* Number of dimensions of the netCDF file */
  int num_dimensions { 0 };         //!< The total number of dimensions defined within the netCDF file
  int num_variables { 0 };          //!< The total number of variables defined within the netCDF file
  int num_global_attributes { 0 };  //!< The total number of global attributes given within the netCDF file
  int id_unlimited_dim { -1 };      //!< The dimension ID of an unlimited dimension (mostly 'Time'-coordinate)

  /* Saves corresponding data to each inquired data variable */
  std::vector<t8_geo_var_t> vars;  //!< A vector holding all variables which has been read from the netCDF file

  t8_nc_mesh_t mesh { nullptr };  //!< A variable holding the constructed forest corresponding to the netCDF data
  //int dimensionality{-1}; //!< The actual dimensionality of the data

  /* These variables save all dimension names and sizes */
  std::vector<size_t> dimension_sizes {};
  std::vector<std::string> dimension_names;

  /* A status flag of the current mode the struct is in */
  bool use_distributed_data {
    false
  };  //!< A flag whether or not the data is/will be distributed among several processes
  sc_MPI_Comm comm { sc_MPI_COMM_WORLD };  //!< The communicator to use in a parallel environment

  t8_nc_par_reading_distribution reading_distribution {
    t8_nc_par_reading_distribution::T8_NC_PAR_DISTRIBUTION_UNDEFINED
  };  //!< Defines the way the netCDF data is read in parallel
  std::vector<int> reading_distribution_num_procs_per_dims;

  int
  get_ncid () const;  //!< Return the id of the corresponding netCDF file
};

int
t8_nc_data::get_ncid () const
{
  return this->ncid;
}

static void
t8_nc_set_dimension_hint_internal (t8_nc_data_t nc_data, std::basic_string<char>&& dim_name,
                                   const t8_coord_ids coord_id)
{
#ifdef T8_WITH_NETCDF
  /* Save the name as a hint for the corresponding dimension */
  nc_data->dimension_hints[dim_name] = coord_id;
#endif
}

void
t8_nc_set_hint_interpret_as_x_axis (t8_nc_data_t nc_data, const char* dimension_name)
{
#ifdef T8_WITH_NETCDF
  /* Set the hint for the longitude dimension */
  t8_nc_set_dimension_hint_internal (nc_data, std::basic_string<char> (dimension_name), t8_coord_ids::T8_LON);
#endif
}

void
t8_nc_set_hint_interpret_as_y_axis (t8_nc_data_t nc_data, const char* dimension_name)
{
#ifdef T8_WITH_NETCDF
  /* Set the hint for the latitude dimension */
  t8_nc_set_dimension_hint_internal (nc_data, std::basic_string<char> (dimension_name), t8_coord_ids::T8_LAT);
#endif
}

void
t8_nc_set_hint_interpret_as_z_axis (t8_nc_data_t nc_data, const char* dimension_name)
{
#ifdef T8_WITH_NETCDF
  /* Set the hint for the vertical dimension */
  t8_nc_set_dimension_hint_internal (nc_data, std::basic_string<char> (dimension_name), t8_coord_ids::T8_LEV);
#endif
}

void
t8_nc_set_hint_interpret_as_time_axis (t8_nc_data_t nc_data, const char* dimension_name)
{
#ifdef T8_WITH_NETCDF
  /* Set the hint for the time dimension */
  t8_nc_set_dimension_hint_internal (nc_data, std::basic_string<char> (dimension_name), t8_coord_ids::T8_TIME);
#endif
}

void
t8_nc_set_hint_read_data_blockwise_in_parallel (t8_nc_data_t nc_data, const int num_dimensions,
                                                int* num_processes_per_dimension)
{
#ifdef T8_WITH_NETCDF_PAR
  std::vector<int> num_procs_per_dim;
  num_procs_per_dim.reserve (num_dimensions);

  /* Transform the array pointer to an vector holding the data */
  for (int i = 0; i < num_dimensions; ++i) {
    num_procs_per_dim.push_back (num_processes_per_dimension[i]);
  }

  /* Store the hint */
  std::swap (nc_data->reading_distribution_num_procs_per_dims, num_procs_per_dim);
  nc_data->reading_distribution = t8_nc_par_reading_distribution::T8_NC_PAR_BLOCKED;
#endif
}

/**
 * \brief Inquire some basic properties of the netCDF file (like the number of variables etc.) 
 * 
 * \param nc_data The \struct t8_nc_data holding the handle to the netCDF file
 */
static void
t8_nc_inquire_general_information (t8_nc_data_t nc_data)
{
#ifdef T8_WITH_NETCDF
  /* Inquire the most general information about the supplied netCDF file (like number of diemnsion, number of global attributes, number of variables, etc.)*/
  int err = nc_inq (nc_data->get_ncid (), &(nc_data->num_dimensions), &(nc_data->num_variables),
                    &(nc_data->num_global_attributes), &(nc_data->id_unlimited_dim));
  t8_nc_check_err (err);
#endif
}

static t8_nc_data_ordering
t8_nc_geo_variable_get_ordering_from_linear_axis_ordering (const std::vector<int>& axis_ordering)
{
#ifdef T8_WITH_NETCDF
  /* Only 2D and 3D geo-spatial data variables are supported for the data layout */
  T8_ASSERT (axis_ordering.size () == 2 || axis_ordering.size () == 3);

  /* Declare a variable holding the computed layout/ordering */
  t8_nc_data_ordering data_layout;

  /* Check the dimensionality of the axis ordering */
  if (axis_ordering.size () == 2) {
    switch (axis_ordering[0]) {
    case t8_coord_ids::T8_LON:
      if (axis_ordering[1] == t8_coord_ids::T8_LAT) {
        data_layout = t8_nc_data_ordering::T8_2D_LON_LAT;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_2D_LON_LEV;
      }
      break;
    case t8_coord_ids::T8_LAT:
      if (axis_ordering[1] == t8_coord_ids::T8_LON) {
        data_layout = t8_nc_data_ordering::T8_2D_LAT_LON;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_2D_LAT_LEV;
      }
      break;
    case t8_coord_ids::T8_LEV:
      if (axis_ordering[1] == t8_coord_ids::T8_LAT) {
        data_layout = t8_nc_data_ordering::T8_2D_LEV_LAT;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_2D_LEV_LON;
      }
      break;
    default:
      t8_errorf ("No geo-spatial coordinate was found in the axis ordering's first position.\n");
    }
  }
  else if (axis_ordering.size () == 3) {
    switch (axis_ordering[0]) {
    case t8_coord_ids::T8_LON:
      if (axis_ordering[1] == t8_coord_ids::T8_LAT) {
        data_layout = t8_nc_data_ordering::T8_3D_LON_LAT_LEV;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_3D_LON_LEV_LAT;
      }
      break;
    case t8_coord_ids::T8_LAT:
      if (axis_ordering[1] == t8_coord_ids::T8_LON) {
        data_layout = t8_nc_data_ordering::T8_3D_LAT_LON_LEV;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_3D_LAT_LEV_LON;
      }
      break;
    case t8_coord_ids::T8_LEV:
      if (axis_ordering[1] == t8_coord_ids::T8_LON) {
        data_layout = t8_nc_data_ordering::T8_3D_LEV_LON_LAT;
      }
      else {
        data_layout = t8_nc_data_ordering::T8_3D_LEV_LAT_LON;
      }
      break;
    default:
      t8_errorf ("No geo-spatial coordinate was found in the axis ordering's first position.\n");
    }
  }

  /* Return the computed layout/ordering of the data */
  return data_layout;
#endif
}

/**
 * \brief Inquire the supplied dimensions
 * 
 * \param nc_data 
 * \param dimensionality 
 * \param dimension_names 
 */
static void
t8_nc_inquire_dimensions (t8_nc_data_t nc_data)
{
#ifdef T8_WITH_NETCDF

  /* Declare a variable for the obtained dimension name */
  std::basic_string<char> nc_name;
  nc_name.reserve (NC_MAX_NAME + 1);

  /* Declare a variable for the obtained dimension length */
  size_t dim_length;

  /* Allocate memory for all 'dimension_lengths' and 'dimension_names' */
  nc_data->dimension_sizes.reserve (nc_data->num_dimensions);
  nc_data->dimension_names.reserve (nc_data->num_dimensions);

  /* Loop over all dimension within the netCDF file */
  for (int dim_index { 0 }; dim_index < nc_data->num_dimensions; ++dim_index) {
    /* Get the name and the length of the dimension which corresponds to the current id */
    int err = nc_inq_dim (nc_data->get_ncid (), dim_index, &nc_name[0], &dim_length);
    t8_nc_check_err (err);

    /* Search for the dimension name within the hints */
    auto dim_hint = nc_data->dimension_hints.find (nc_name.c_str ());

    /* Check if a hint was found */
    if (dim_hint != nc_data->dimension_hints.end ()) {
      /* Save the id of this dimension at the position of the hint */
      nc_data->coord_dim_ids[dim_hint->second] = dim_index;

      /* Save the global length of the dimension */
      nc_data->coord_lengths[dim_hint->second] = dim_length;
    }
    else {
      /* In case no hint was found corresponding to this name, we will check it for some standard names */
      if (nc_data->coord_dim_ids[t8_coord_ids::T8_LON] == T8_NC_NOT_CONSIDERED
          && (strcmp (nc_name.c_str (), "lon") == 0 || strcmp (nc_name.c_str (), "longitude") == 0)) {
        /* Save the netCDF internal ID of this diemnsion and the global length of it */
        nc_data->coord_dim_ids[t8_coord_ids::T8_LON] = dim_index;
        nc_data->coord_lengths[t8_coord_ids::T8_LON] = dim_length;
      }
      else if (nc_data->coord_dim_ids[t8_coord_ids::T8_LAT] == T8_NC_NOT_CONSIDERED
               && (strcmp (nc_name.c_str (), "lat") == 0 || strcmp (nc_name.c_str (), "latitude") == 0)) {
        /* Save the netCDF internal ID of this diemnsion and the global length of it */
        nc_data->coord_dim_ids[t8_coord_ids::T8_LAT] = dim_index;
        nc_data->coord_lengths[t8_coord_ids::T8_LAT] = dim_length;
      }
      else if (nc_data->coord_dim_ids[t8_coord_ids::T8_LEV] == T8_NC_NOT_CONSIDERED
               && (strcmp (nc_name.c_str (), "lev") == 0 || strcmp (nc_name.c_str (), "height") == 0
                   || strcmp (nc_name.c_str (), "depth") == 0)) {
        /* Save the netCDF internal ID of this diemnsion and the global length of it */
        nc_data->coord_dim_ids[t8_coord_ids::T8_LEV] = dim_index;
        nc_data->coord_lengths[t8_coord_ids::T8_LEV] = dim_length;
      }
      else if (nc_data->coord_dim_ids[t8_coord_ids::T8_TIME] == T8_NC_NOT_CONSIDERED
               && strcmp (nc_name.c_str (), "time") == 0) {
        /* Save the netCDF internal ID of this diemnsion and the global length of it */
        nc_data->coord_dim_ids[t8_coord_ids::T8_TIME] = dim_index;
        nc_data->coord_lengths[t8_coord_ids::T8_TIME] = dim_length;
      }
    }

    /* Save the name and the length of each dimension within the netCDF file */
    nc_data->dimension_names.push_back (std::string (nc_name.c_str ()));
    nc_data->dimension_sizes.push_back (dim_length);
  }

#endif
}

/* Inquire the coordinate dimension data */
static void
t8_nc_inquire_dimensions_data (t8_nc_data_t nc_data)
{
#ifdef T8_WITH_NETCDF
  int err;
  nc_type var_type;

  /* Loop over all possibly-considered coordinate variables (latitude, longitude, height, time) */
  for (int coord_id { 0 }; coord_id < t8_coord_ids::T8_NUM_COORD_IDS; ++coord_id) {
    /* Check if the coordinate variable/dimension is supplied/considered within the netCDF file */
    if (nc_data->coord_dim_ids[coord_id] != T8_NC_NOT_CONSIDERED) {
      /* Coordinate variables have a concerning coordinate dimension; the name of the dimension coincides with the name of the variable and is only dependent on "its own dimension" */
      /* Get the ID of the coordinate variable */
      err = nc_inq_varid (nc_data->get_ncid (), nc_data->dimension_names[nc_data->coord_dim_ids[coord_id]].c_str (),
                          &(nc_data->coord_var_ids[coord_id]));
      t8_nc_check_err (err);

      /* Get the type of coordinate variables (e.g. float or int) */
      err = nc_inq_vartype (nc_data->get_ncid (), nc_data->coord_var_ids[coord_id], &(var_type));
      t8_nc_check_err (err);

      /* Allocate space for each (geo-spatial) coordinate variable */
      nc_data->coordinates[coord_id] = t8_nc_create_geo_variable (
        nc_data->dimension_names[nc_data->coord_dim_ids[coord_id]].c_str (), static_cast<t8_geo_data_type> (var_type),
        nc_data->dimension_sizes[nc_data->coord_dim_ids[coord_id]]);

      /* Read in the data from the coordinate variable */
      err = nc_get_var (nc_data->get_ncid (), nc_data->coord_var_ids[coord_id],
                        t8_nc_geo_variable_get_data_ptr (nc_data->coordinates[coord_id]));
      t8_nc_check_err (err);

      t8_global_productionf ("Coordinate data of variable %s has been inquired.\n",
                             nc_data->dimension_names[nc_data->coord_dim_ids[coord_id]].c_str ());
    }
  }

#endif
}

/**
 * \brief A struct for handling input data from netCDF files. This struct holds
 * netCDF variable related information only
 */
struct t8_netcdf_var_data
{
 private:
  int var_id { -1 };  //!< The netCDF file internal id of this variable
 public:
  t8_netcdf_var_data (const int _var_id, const int _var_type): var_id { _var_id }, var_type { _var_type } {};

  int var_type { -1 };             //!< NetCDF variable type
  std::vector<int> dimension_ids;  //!< The ids of the dimension this variable consists of
  /* Describing the netCDF data hyperslab */
  std::vector<size_t> start_ptr;  //!< Coordinate start values of the data
  std::vector<size_t> count_ptr;  //!< Coordinate length values of the data

  /* Axis-ordering of the data */
  std::vector<int> axis_ordering;  //!< Specifies the axis-ordering of the variable

  /* Function for getting the netCDF variable ID */
  int
  get_var_id () const
  {
    return var_id;
  };
};

/* Typedef for a pointer to a struct */
typedef struct t8_netcdf_var_data* t8_netcdf_var_data_t;

static void
t8_nc_constrain_dimensions_to_given_hyperslab (t8_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr)
{
#ifdef T8_WITH_NETCDF
  /* In order to restrict the dimensions to a hyperslab, we need at least one variable to which this hyperslab corresponds to */
  T8_ASSERT (nc_data->vars.size () > 0);

  /* Get the netCDF specific variable data of the first variable */
  t8_netcdf_var_data_t nc_var_data
    = static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (nc_data->vars.front ()));

  /* Eventually adjust the global coordinates */
  for (size_t dims { 0 }; dims < nc_var_data->dimension_ids.size (); ++dims) {
    /* Get the dimension id of this variable */
    const int current_dim_id = nc_var_data->dimension_ids[dims];

    /* Check to which coordinate dimension the dimension ID complies */
    if (current_dim_id == nc_data->coord_dim_ids[t8_coord_ids::T8_LON]) {
      if (count_ptr[dims] < t8_nc_geo_variable_get_num_elements (nc_data->coordinates[t8_coord_ids::T8_LON])) {
        /* Only a part of the coordinate dimension is considered */
        t8_nc_geo_variable_crop_data_to_selection (nc_data->coordinates[t8_coord_ids::T8_LON], start_ptr[dims],
                                                   start_ptr[dims] + count_ptr[dims] - 1);
        /* Set the coordinate length correspondingly */
        nc_data->coord_lengths[t8_coord_ids::T8_LON] = count_ptr[dims];
      }
    }
    else if (current_dim_id == nc_data->coord_dim_ids[t8_coord_ids::T8_LAT]) {
      if (count_ptr[dims] < t8_nc_geo_variable_get_num_elements (nc_data->coordinates[t8_coord_ids::T8_LAT])) {
        /* Only a part of the coordinate dimension is considered */
        t8_nc_geo_variable_crop_data_to_selection (nc_data->coordinates[t8_coord_ids::T8_LAT], start_ptr[dims],
                                                   start_ptr[dims] + count_ptr[dims] - 1);
        /* Set the coordinate length correspondingly */
        nc_data->coord_lengths[t8_coord_ids::T8_LAT] = count_ptr[dims];
      }
    }
    else if (current_dim_id == nc_data->coord_dim_ids[t8_coord_ids::T8_LEV]) {
      if (count_ptr[dims] < t8_nc_geo_variable_get_num_elements (nc_data->coordinates[t8_coord_ids::T8_LEV])) {
        /* Only a part of the coordinate dimension is considered */
        t8_nc_geo_variable_crop_data_to_selection (nc_data->coordinates[t8_coord_ids::T8_LEV], start_ptr[dims],
                                                   start_ptr[dims] + count_ptr[dims] - 1);
        /* Set the coordinate length correspondingly */
        nc_data->coord_lengths[t8_coord_ids::T8_LEV] = count_ptr[dims];
      }
    }
    else if (current_dim_id == nc_data->coord_dim_ids[t8_coord_ids::T8_TIME]) {
      if (count_ptr[dims] < t8_nc_geo_variable_get_num_elements (nc_data->coordinates[t8_coord_ids::T8_TIME])) {
        /* Only a part of the coordinate dimension is considered */
        t8_nc_geo_variable_crop_data_to_selection (nc_data->coordinates[t8_coord_ids::T8_TIME], start_ptr[dims],
                                                   start_ptr[dims] + count_ptr[dims] - 1);
        /* Set the coordinate length correspondingly */
        nc_data->coord_lengths[t8_coord_ids::T8_TIME] = count_ptr[dims];
      }
    }
  }
#endif
}

static void
t8_nc_constrain_dimensions_to_given_dimensionality (t8_nc_data_t nc_data, const int dimensionality)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (dimensionality > 0 && dimensionality <= 4);
  /* Iterate over all geo-spatial dimensions and check whether they are supplied */
  for (int i = 0; i < static_cast<int> (t8_coord_ids::T8_NUM_COORD_IDS); ++i) {
    if (i < dimensionality) {
      /* Check whether an id has been found */
      if (nc_data->coord_dim_ids[i] == T8_NC_NOT_CONSIDERED) {
        /* If no ID is associated with this dimension, we cannot built the mesh based on the supplied dimensionality */
        t8_errorf ("The mesh cannot be built based on the given dimensionality (= %d), because there was no ... "
                   "dimension found.\nBuilding the mesh requires a longitude, latitude and vertical coordinate in the "
                   "given order for 1D, 2D or 3D meshes. If another order is wished for or different dimension names "
                   "are present, please use the t8_nc_set_hint_interpret_as_..._axis-functions.\n",
                   dimensionality);
      }
    }
    else {
      /* All other dimensions which were found but will not be considered, will be reset */
      nc_data->coord_dim_ids[i] = T8_NC_NOT_CONSIDERED;
      nc_data->coord_lengths[i] = 0;
    }
  }
#endif
}

/** Check if given attributes are stored alongside the inquired variables
 * \note Only if the data type of the attribute is (signed/unsigned) short, int, int64, double or float, the data value will be stored (and applied later on) */
static void
t8_nc_check_var_for_attributes (t8_nc_data_t nc_data, t8_geo_var_t var)
{
#ifdef T8_WITH_NETCDF
  int err;
  int attribute_type;
  t8_universal_type_t att_value;

  /* The attributes which will be checked */
  const std::vector<const char*> att_names { "add_offset", "scale_factor", "missing_value" };

  /* Obtain the netCDF specific variable data stored within the user data */
  t8_netcdf_var_data_t nc_var_data = static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (var));

  t8_global_productionf ("Checking attributes of variable %s...\n", t8_nc_geo_variable_get_name (var));

  /* Iterate over all attributes */
  for (auto iter_att { att_names.begin () }; iter_att != att_names.end (); ++iter_att) {
    /* Check the data type of the attribute */
    err = nc_inq_atttype (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &attribute_type);

    if (err == NC_NOERR) {
      t8_global_productionf ("Found a/an '%s' attribute.\n", *iter_att);

      /* If the data type of the attribute has been successfully inquired, get the data value of the attribute */
      switch (attribute_type) {
      case NC_SHORT:
        int16_t short_att;
        err = nc_get_att_short (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &short_att);
        t8_nc_check_err (err);
        att_value = static_cast<int16_t> (short_att);
        break;
      case NC_INT:
        int32_t int_att;
        err = nc_get_att_int (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &int_att);
        t8_nc_check_err (err);
        att_value = static_cast<int32_t> (int_att);
        break;
      case NC_FLOAT:
        float float_att;
        err = nc_get_att_float (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &float_att);
        t8_nc_check_err (err);
        att_value = static_cast<float> (float_att);
        break;
      case NC_DOUBLE:
        double double_att;
        err = nc_get_att_double (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &double_att);
        t8_nc_check_err (err);
        att_value = static_cast<double> (double_att);
        break;
      case NC_USHORT:
        uint16_t ushort_att;
        err = nc_get_att_ushort (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &ushort_att);
        t8_nc_check_err (err);
        att_value = static_cast<uint16_t> (ushort_att);
        break;
      case NC_UINT:
        uint32_t uint_att;
        err = nc_get_att_uint (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &uint_att);
        t8_nc_check_err (err);
        att_value = static_cast<uint32_t> (uint_att);
        break;
      case NC_INT64:
        long long int longlong_att;
        err = nc_get_att_longlong (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &longlong_att);
        t8_nc_check_err (err);
        att_value = static_cast<int64_t> (longlong_att);
        break;
      case NC_UINT64:
        unsigned long long int ulonglong_att;
        err = nc_get_att_ulonglong (nc_data->get_ncid (), nc_var_data->get_var_id (), *iter_att, &ulonglong_att);
        t8_nc_check_err (err);
        att_value = static_cast<uint64_t> (ulonglong_att);
        break;
      default:
        t8_global_productionf ("The attribute  %s of variable %s could not be inquired (and will not be applied), "
                               "since its data type is not supported.\n",
                               *iter_att, t8_nc_geo_variable_get_name (var));
      }

      /* Check which possible attribute is considered and save it's value */
      if (strcmp (*iter_att, "add_offset") == 0) {
        t8_nc_geo_variable_set_add_offset (var, att_value);
      }
      else if (strcmp (*iter_att, "scale_factor") == 0) {
        t8_nc_geo_variable_set_scale_factor (var, att_value);
      }
      else if (strcmp (*iter_att, "missing_value") == 0) {
        t8_nc_geo_variable_set_missing_value (var, att_value);
      }
    }
  }
#endif
};

static void
t8_nc_inquire_variables (t8_nc_data_t nc_data, const int num_variables, const char** const variable_names)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (num_variables > 0);

  int err, nc_internal_var_id, num_dimensions;
  nc_type var_type;
  std::array<int, NC_MAX_VAR_DIMS> dim_ids;

  std::vector<t8_netcdf_var_data> nc_var_info;
  nc_var_info.reserve (num_variables);

  /* Iterate over all variables */
  for (int var_id = 0; var_id < num_variables; ++var_id) {
    /* Inquire the variable id of the given variable's name */
    err = nc_inq_varid (nc_data->get_ncid (), variable_names[var_id], &(nc_internal_var_id));
    if (err == NC_NOERR) {
      /* If an ID is found, information like number of dimensions, data type, etc. will be inquired */
      err = nc_inq_var (nc_data->get_ncid (), nc_internal_var_id, NULL, &(var_type), &(num_dimensions), dim_ids.data (),
                        NULL);
      t8_nc_check_err (err);

      /* Create a new variable with the name and the data type */
      nc_data->vars.push_back (
        t8_nc_create_geo_variable (variable_names[var_id], t8_nc_geo_variable_nc_type_to_t8_geo_data_type (var_type)));

      /* Create a new netCDF specific variable data struct */
      t8_netcdf_var_data_t nc_var_data = new t8_netcdf_var_data (nc_internal_var_id, var_type);

      /* Allocate memory for the dimension ids this variable depends on */
      nc_var_data->dimension_ids.reserve (num_dimensions);

      /* Store the ids within the netCDF specific variable struct */
      std::copy_n (dim_ids.begin (), num_dimensions, std::back_inserter (nc_var_data->dimension_ids));

      /* Assign the netCDF specific variable data to the newly created variable */
      t8_nc_geo_variable_set_user_data (nc_data->vars.back (), static_cast<void*> (nc_var_data));

      /* Check if an 'offset'-attribute, 'missing_value'-attribute or 'scale_factor'-attribute is given for each variable */
      t8_nc_check_var_for_attributes (nc_data, nc_data->vars.back ());
    }
    else {
      /* Print an error message if no ID to the supplied name of the variable is found */
      t8_global_productionf ("An error occurred while inquiring the netCDf 'variable_id' corresponding to the supplied "
                             "'variable_name' (%s).\n",
                             variable_names[var_id]);
      t8_nc_check_err (err);
    }
  }
#endif
}

/**
 * \brief Open a netCDF file and return a struct holding the handle to it 
 */
t8_nc_data_t
t8_nc_start (const char* path_to_file, const enum t8_nc_opening_mode mode, const sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
  int err { 0 };
  int ncid;
  bool parallel_access { false };

  switch (mode) {
  case t8_nc_opening_mode::T8_NC_SERIAL:
    /* Open the file in a serial mode */
    ncid = t8_nc_open_serial (path_to_file);
    break;
  case t8_nc_opening_mode::T8_NC_PARALLEL: {
#ifdef T8_WITH_NETCDF_PAR
    int comm_size { 0 };
    /* Get the size of the supplied communicator */
    err = sc_MPI_Comm_size (comm, &comm_size);
    SC_CHECK_MPI (err);

    /* Check whether there are several processes in the communicator for actually reading in parallel from the file */
    if (comm_size > 1) {
      /* Open the file for parallel access */
      ncid = t8_nc_open_parallel (path_to_file, comm);

      /* Set the flag for parallel access */
      parallel_access = true;
    }
    else {
      t8_debugf ("The file is ought to be opened in a parallel mode but there is not more than one process within the "
                 "supplied communicator. Therefore, the file will be opened in a serial mode.\n");

      /* Open the file in a serial mode */
      ncid = t8_nc_open_serial (path_to_file);
    }
#else
    t8_errorf ("NetCDF's parallel file access functions are not available. Please ensure that 'netcdf_par.h' is "
               "present when t8code is linked against netCDF.");
#endif
  } break;
  default:
    t8_errorf ("An unknown netCDF opening mode was supplied.");
  }

  /* Create a netCDF data struct and return it */
  t8_nc_data_t nc_data = new t8_nc_data { ncid, parallel_access, comm };

  /* Inquire some general information about the file */
  t8_nc_inquire_general_information (nc_data);

  return nc_data;

#else
  t8_errorf ("t8code is not linked against netCDF. Please recompile the build with netCDF linkage.\n");
  return T8_NC_ERR_CODE;
#endif
}

/**
 * \brief Opens a file explicitly for serial read access
 */
int
t8_nc_open_serial (const char* path_to_file)
{
#ifdef T8_WITH_NETCDF
  int err { 0 };
  int ncid { 0 };

  /* Open the file without explicit parallel access */
  err = nc__open (path_to_file, NC_NOWRITE, NULL, &ncid);
  t8_nc_check_err (err);

  t8_global_productionf ("The netCDF-File %s has been opened for serial reading.\n", path_to_file);

  return ncid;
#else
  t8_errorf ("t8code is not linked against netCDF. Please recompile the build with netCDF\n.");
  return T8_NC_ERR_CODE;
#endif
}

/**
 * \brief Opens a file explicitly for parallel read access 
 */
int
t8_nc_open_parallel (const char* path_to_file, MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF_PAR
  int err { 0 };
  int ncid { 0 };

  int comm_size { 0 };
  err = sc_MPI_Comm_size (comm, &comm_size);
  SC_CHECK_MPI (err);

  /* Is more than one process present */
  T8_ASSERT (comm_size > 1);

  sc_MPI_Info info = sc_MPI_INFO_NULL;
  err = nc_open_par (path_to_file, NC_NOWRITE, comm, info, &ncid);
  t8_nc_check_err (err);

  t8_global_productionf ("The netCDF-File %s has been opened for parallel reading.\n", path_to_file);

  return ncid;
#else
  t8_errorf ("t8code is not linked against netCDF (at least not to it's parallel functionality). Please recompile the "
             "build with netCDF and ensure accessibility to the 'netcdf_par.h' header file.\n.");
  return T8_NC_ERR_CODE;
#endif
}

/**
 * \brief Close the netCDF file
 */
void
t8_nc_close (int ncid)
{
#ifdef T8_WITH_NETCDF
  int err { nc_close (ncid) };
  t8_nc_check_err (err);

  t8_global_productionf ("The netCDF-File has been closed.\n");
#endif
}

/**
 * \brief This function closes the netCDF file and deallocates the internal data structure
 */
void
t8_nc_finish (t8_nc_data_t nc_data)
{
#ifdef T8_WITH_NETCDF
  /* Close the netCDF File */
  t8_nc_close (nc_data->get_ncid ());

  /* Deallocate the nc_data */
  if (nc_data != nullptr) {
    delete nc_data;
  }
#else
  t8_errorf ("t8code is not linked against netCDF. Please recompile the build with netCDF\n.");
#endif
}

/* Calculate a blocked distribution based on supplied numbers for the amount of processes per dimension */
static std::pair<std::vector<size_t>, std::vector<size_t>>
t8_nc_calculate_reading_data_distribution_num_procs_per_dim (const std::vector<size_t>& start_values,
                                                             const std::vector<size_t>& count_values,
                                                             const sc_MPI_Comm comm,
                                                             const std::vector<int>& num_procs_per_dim)
{
  T8_ASSERT (start_values.size () == count_values.size () && count_values.size () == num_procs_per_dim.size ());
  int err, rank, size;

  /* Get the size of the communicator */
  err = sc_MPI_Comm_size (comm, &size);
  SC_CHECK_MPI (err);

  /* Check if the blocked distribution given by the parameter @var num_procs_per_dim equals or is less than the size of the communicator */
  T8_ASSERT (std::accumulate (num_procs_per_dim.begin (), num_procs_per_dim.end (), 1, std::multiplies<int> ())
             <= size);

  /* Get the (local) rank id */
  err = sc_MPI_Comm_rank (comm, &rank);
  SC_CHECK_MPI (err);

  /* Get the amount of processes reading from the netCDF file */
  const int num_reading_procs
    = std::accumulate (num_procs_per_dim.begin (), num_procs_per_dim.end (), 1, std::multiplies<int> ());

  /* Allocate vectors holding the local start and count offsets */
  std::vector<size_t> start_vals;
  start_vals.reserve (start_values.size ());
  std::vector<size_t> count_vals;
  count_vals.reserve (start_values.size ());

  /* Check if the process-local rank will read data */
  if (rank < num_reading_procs) {
    /* Define a variable keeping track of the current index in the vectors */
    int current_dim_id { 0 };

    /* Variables for calculating the distribution */
    size_t global_offset { 0 };
    size_t box_length { 0 };

    /* Variabels for calculating a 'box'-ed hyperslab*/
    int vrank = rank;
    int divisor = 1;
    size_t div = 1;

    /* Iterate over all dimensions and determine a blocked distribution */
    for (auto iter { count_values.begin () }; iter != count_values.end (); ++iter, ++current_dim_id) {
      T8_ASSERT (num_procs_per_dim[current_dim_id] <= num_reading_procs);

      /* If a dimension is not considered in the data hyperslab, just skip it */
      if (*iter <= 1UL) {
        /* These start and count values indicate that we are not reading data from this dimension */
        start_vals.push_back (start_values[current_dim_id]);
        count_vals.push_back (1);
      }
      else {
        /* Check if there are more than one rank considered for this dimension */
        if (num_procs_per_dim[current_dim_id] > 1) {
          if (div > 1) {
            vrank = rank / std::pow (2, div - 1);
          }
          else {
            vrank = rank;
          }
          divisor = num_procs_per_dim[current_dim_id];
          global_offset
            = static_cast<uint64_t> ((((double) (vrank % divisor) * (long double) (*iter)) / (double) divisor));

          box_length = (vrank % divisor) + 1 >= num_procs_per_dim[current_dim_id]
                         ? (*iter) - global_offset
                         : static_cast<uint64_t> (
                             (((double) ((vrank % divisor) + 1) * (long double) (*iter)) / (double) divisor))
                             - global_offset;

          ++div;
        }
        else {
          global_offset = 0;
          box_length = *iter;
        }

        /* Save the process-local start and count values for the current dimension */
        start_vals.push_back (start_values[current_dim_id] + global_offset);
        count_vals.push_back (box_length);
      }
    }
  }
  else {
    for (auto iter { count_values.begin () }; iter != count_values.end (); ++iter) {
      start_vals.push_back (0UL);
      count_vals.push_back (0UL);
    }
  }

  /* Return a struct dscribing the data which should be read in on each rank */
  return std::make_pair (std::move (start_vals), std::move (count_vals));
}

std::pair<std::vector<size_t>, std::vector<size_t>>
t8_nc_calculate_reading_data_distribution_blocked (const std::vector<size_t>& start_values,
                                                   const std::vector<size_t>& count_values, const MPI_Comm comm,
                                                   const std::vector<int>& num_procs_per_dim)
{
  return t8_nc_calculate_reading_data_distribution_num_procs_per_dim (start_values, count_values, comm,
                                                                      num_procs_per_dim);
}

/** Inquire the data for each given variable (either the variable as a whole (start_ptr and count_ptr equal to nullptrs) or a specified hyperslab) */
static void
t8_nc_inquire_variables_data (t8_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr)
{
#ifdef T8_WITH_NETCDF

  int err;
  size_t num_data_points { 1 };

  /* Vectors storing the start and count values for reading the data */
  std::vector<size_t> start_values;
  std::vector<size_t> count_values;

  /** Iterate over all variables and save their axis-ordering 
      * (since the specified hyperslab corresponds to all given variables, their axis-ordering should coincide) */
  for (auto iter { nc_data->vars.begin () }; iter != nc_data->vars.end (); ++iter) {
    /* Get the netCDF specific variable data */
    t8_netcdf_var_data_t nc_var_data = static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (*iter));

    /** Since all variables are defined by the same hyperslab, we can assume that the have the same axis-ordering 
         * Therefore, we are computing the (reduced) axis-ordering only for the first variable and copy the ordering for all remaining variables.
         */

    /* Get the iterator to the first variable */
    auto iter_to_first_variable = nc_data->vars.begin ();

    /* Check if we are in the first iteration. In this case the axis ordering needs to be calculated otherwise just copied */
    if (iter == iter_to_first_variable) {
      /* Define a new axis-ordering vector */
      std::vector<int> axis_ordering;
      axis_ordering.reserve (t8_coord_ids::T8_NUM_COORD_IDS);

      /* Define a variable for accessing the hyperslab specification arrays */
      int hyperslab_access_id = 0;

      /* Retrieve the axis ordering of the variable */
      for (auto dim_iter { nc_var_data->dimension_ids.begin () }; dim_iter != nc_var_data->dimension_ids.end ();
           ++dim_iter, ++hyperslab_access_id) {
        if (*dim_iter == nc_data->coord_dim_ids[t8_coord_ids::T8_LON]) {
          if (count_ptr[hyperslab_access_id] > 1) {
            axis_ordering.push_back (t8_coord_ids::T8_LON);
          }
        }
        else if (*dim_iter == nc_data->coord_dim_ids[t8_coord_ids::T8_LAT]) {
          if (count_ptr[hyperslab_access_id] > 1) {
            axis_ordering.push_back (t8_coord_ids::T8_LAT);
          }
        }
        else if (*dim_iter == nc_data->coord_dim_ids[t8_coord_ids::T8_LEV]) {
          if (count_ptr[hyperslab_access_id] > 1) {
            axis_ordering.push_back (t8_coord_ids::T8_LEV);
          }
        }
        else if (*dim_iter == nc_data->coord_dim_ids[t8_coord_ids::T8_TIME]) {
          //Time series are currently not supported
          if (count_ptr[hyperslab_access_id] > 1) {
            t8_errorf ("Time-Series data is not yet supported, please choose a single point in time.\n");
            axis_ordering.push_back (t8_coord_ids::T8_TIME);
          }
        }

        /* We are 'transforming' the hyperslab pointers to vectors */
        start_values.push_back (start_ptr[hyperslab_access_id]);
        count_values.push_back (count_ptr[hyperslab_access_id]);
      }

      /* Save the axis ordering */
      nc_var_data->axis_ordering = std::move (axis_ordering);

      /* Get the data ordering (enumerator) corresponding to this axis ordering */
      t8_nc_geo_variable_set_data_ordering_scheme (
        *iter, t8_nc_geo_variable_get_ordering_from_linear_axis_ordering (nc_var_data->axis_ordering));
    }
    else {
      /* Get the netCDF specific variable data of the first variable (which has already stored the axis ordering) */
      t8_netcdf_var_data_t nc_first_var_data
        = static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (*iter_to_first_variable));

      /* Copy the axis ordering */
      nc_var_data->axis_ordering = nc_first_var_data->axis_ordering;

      /* Copy the data ordering */
      t8_nc_geo_variable_set_data_ordering_scheme (
        *iter, t8_nc_geo_variable_get_data_ordering_scheme (*iter_to_first_variable));
    }

    /* Check if we are in the first iteration. In this case the (local) hyperslab of a variable needs to be calculated otherwise just copied */
    if (iter == iter_to_first_variable) {
      /* Check if the variables' data may be read in distributed */
      if (nc_data->use_distributed_data) {
        /* Declare a pair which will hold the computed start and count offsets */
        std::pair<std::vector<size_t>, std::vector<size_t>> local_hyperslab;

        /* Calculate the offsets for each process based on the (preferred) data distribution */
        switch (nc_data->reading_distribution) {
        case t8_nc_par_reading_distribution::T8_NC_PAR_BLOCKED:
          /* Calculate a blocked data dsitribution */
          local_hyperslab = t8_nc_calculate_reading_data_distribution_blocked (
            start_values, count_values, nc_data->comm, nc_data->reading_distribution_num_procs_per_dims);
          break;
        case t8_nc_par_reading_distribution::T8_NC_PAR_DISTRIBUTION_UNDEFINED: {
          /* If no distribution has been specified, we default to a blocked distribution */
          int err, comm_size;

          /* Get the size of the MPI communicator */
          err = sc_MPI_Comm_size (nc_data->comm, &comm_size);
          SC_CHECK_MPI (err);

          int num_considered_dims = 0;

          /* Define a vector for a blocked distribution (and initialize all entries to one) */
          std::vector<int> p_distribution (start_values.size (), 1);

          /* Iterate over the count vector and inquire the general data dimensions (e.g. 2D or 3D) */
          for (auto citer { count_values.begin () }; citer != count_values.end (); ++citer) {
            if (*citer > 1) {
              ++num_considered_dims;
            }
          }

          T8_ASSERT (num_considered_dims > 0);

          /* Calculate an equal distribution per dimension */
          const int equal_dim_distribution = static_cast<int> (std::pow (comm_size, 1.0 / num_considered_dims));

          /* Iterate again over the count vector in order to assign a blocked distribution at the correct dimensions */
          int dim_index = 0;

          for (auto citer { count_values.begin () }; citer != count_values.end (); ++citer, ++dim_index) {
            if (*citer > 1) {
              p_distribution[dim_index] = equal_dim_distribution;
            }
          }

          /* Calculate a blocked data dsitribution */
          local_hyperslab = t8_nc_calculate_reading_data_distribution_blocked (start_values, count_values,
                                                                               nc_data->comm, p_distribution);

          /* Set that we have defaulted to a blcoekd scheme */
          nc_data->reading_distribution = t8_nc_par_reading_distribution::T8_NC_PAR_BLOCKED;
        } break;
        default:
          t8_errorf ("An unknown data distribution for parallel reading was supplied which is not yet implemented.\n");
        }

        /* Update the start and count values for the variable to their local view */
        std::swap (nc_var_data->start_ptr, local_hyperslab.first);
        std::swap (nc_var_data->count_ptr, local_hyperslab.second);
      }
      else {
        /* Update the start and count values for the variable in case of a serial read) */
        std::swap (nc_var_data->start_ptr, start_values);
        std::swap (nc_var_data->count_ptr, count_values);
      }

      /* Iterate over the (local) count vector in order to calculate the amount of data points for this variable */
      for (auto citer { nc_var_data->count_ptr.begin () }; citer != nc_var_data->count_ptr.end (); ++citer) {
        /* Reduce the data points per dimension in order to obtain the overall amount of data points of the variable */
        num_data_points *= *citer;
      }
    }
    else {
      /* In case we are not concerned with the first variable, we will just copy the calculated (local) hyperslab for each additional variable */
      /* Get the netCDF specific variable data of the first variable (which has already stored the start and count values) */
      t8_netcdf_var_data_t nc_first_var_data
        = static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (*iter_to_first_variable));

      /* Copy the start and count values */
      nc_var_data->start_ptr = nc_first_var_data->start_ptr;
      nc_var_data->count_ptr = nc_var_data->count_ptr;
    }

    /* Set the number of elements of this variable */
    const size_t num_elements = num_data_points;

    /* Allocate memory for an array able to hold the data of the variable */
    t8_nc_geo_variable_allocate_data (*iter, t8_nc_geo_variable_get_type (*iter), num_elements);

    for (size_t j = 0; j < nc_var_data->start_ptr.size (); ++j) {
      t8_debugf ("At pos: %ld, start vec: %ld; count vec: %ld\n", j, nc_var_data->start_ptr[j],
                 nc_var_data->count_ptr[j]);
    }

    /* Read in the (process-local) hyperslab of the data */
    err = nc_get_vara (nc_data->get_ncid (), nc_var_data->get_var_id (), nc_var_data->start_ptr.data (),
                       nc_var_data->count_ptr.data (), t8_nc_geo_variable_get_data_ptr (*iter));
    t8_nc_check_err (err);

    t8_global_productionf ("The data of variable %s has been inquired.\n", t8_nc_geo_variable_get_name (*iter));

    /* Store the axis ordering of the data */
    if (iter == iter_to_first_variable) {
      /* In case of the first variable, the data ordering has to be computed */
    }
    else {
      /* Since all variables are defined on the same hyperslab the data ordering should coincide with the ordering of the first variable and therefore, can just be copied */
    }
  }

#endif
}

static void
t8_nc_deallocate_var_specific_data (t8_nc_data_t nc_data)
{
  for (auto iter { nc_data->vars.begin () }; iter != nc_data->vars.end (); ++iter) {
    if (t8_nc_geo_variable_get_user_data (*iter) != nullptr) {
      delete static_cast<t8_netcdf_var_data*> (t8_nc_geo_variable_get_user_data (*iter));
    }
  }
}

static void
t8_nc_build_initial_mesh (t8_nc_data_t nc_data, const enum t8_nc_geo_mesh_type mesh_type,
                          const enum t8_nc_geo_mesh_form mesh_form, const enum t8_nc_geo_mesh_elements mesh_elems)
{
#ifdef T8_WITH_NETCDF
  /* Variable for the dimensionality of the mesh */
  int count_considered_dims = 0;

  /* Determine the dimensionality of the data */
  for (auto iter { nc_data->coord_lengths.begin () }; iter != nc_data->coord_lengths.end (); ++iter) {
    /* Check if the coordinate dimension is considered at all */
    if (*iter > 1) {
      /* Count all dimensions with length greater than one */
      ++count_considered_dims;
    }
  }

  /* Allocate a new 't8_nc_mesh' struct which will hold the constructed mesh */
  nc_data->mesh = t8_nc_mesh_create ();

  /* Save the dimensionality of the mesh to be constructed */
  t8_nc_mesh_set_dimensionality (nc_data->mesh, count_considered_dims);

  /* Set the data ordering within the 'nc_mesh' struct (this corresponds to the axis ordering of the variables which is equal for all variables) */
  t8_nc_mesh_set_data_ordering_scheme (nc_data->mesh,
                                       t8_nc_geo_variable_get_data_ordering_scheme (nc_data->vars.front ()));

  /* Set the dimension lengths of the data */
  t8_nc_mesh_set_longitude_length (nc_data->mesh, static_cast<int> (nc_data->coord_lengths[t8_coord_ids::T8_LON]));
  t8_nc_mesh_set_latitude_length (nc_data->mesh, static_cast<int> (nc_data->coord_lengths[t8_coord_ids::T8_LAT]));
  t8_nc_mesh_set_vertical_length (nc_data->mesh, static_cast<int> (nc_data->coord_lengths[t8_coord_ids::T8_LEV]));

  /* Declare a forest variable */
  t8_forest_t forest;

  /* Check first which mesh form should be built (Currently, only a rectangular mesh is possible) */
  switch (mesh_form) {
  case t8_nc_geo_mesh_form::T8_NC_RECTANGULAR:
    /* In case a rectangular mesh should be built */
    /* Check the type of the mesh which is about to be built */
    switch (mesh_type) {
    case t8_nc_geo_mesh_type::T8_NC_EMBEDDED_MESH:
      /* In case an embedded mesh has been chosen (embedding the geo-spatial domain into a larger mesh) */
      forest = t8_nc_build_initial_rectangular_embedded_minimal_mesh (nc_data->mesh, nc_data->comm);
      break;
    case t8_nc_geo_mesh_type::T8_NC_CONGRUENT_MESH:
      /* In case a 'congruent' mesh has been chosen (resembling only and fully the geo-spatial domain of the data) */
      forest = t8_nc_build_initial_rectangular_congruent_mesh (nc_data->mesh, nc_data->comm);
      break;
    default:
      t8_errorf ("A not supported mesh type has been selected for the geo-spatial netCDF data. Please, choose either "
                 "an 'embedded' or a 'congruent' mesh type.\n");
    }
    break;
  case t8_nc_geo_mesh_form::T8_NC_SPHERICAL:
    t8_errorf ("Building a spherical mesh is not yet implemented. Please, choose a 'rectangular' mesh.\n");
    break;
  }

#endif
}

/**
 * \brief Construct a mesh based on the supplied coordinate dimension names
 * 
 * \param[in,out] nc_data 
 * \param[in] dimensionality 
 * \param[in] mesh_type 
 * \param[in] mesh_form 
 * \param[in] mesh_elems 
 */
void
t8_nc_construct_mesh (t8_nc_data_t nc_data, const int dimensionality, const enum t8_nc_geo_mesh_type mesh_type,
                      const enum t8_nc_geo_mesh_form mesh_form, const enum t8_nc_geo_mesh_elements mesh_elems)
{
#ifdef T8_WITH_NETCDF

  /* Inquire the coordinate dimensions information */
  t8_nc_inquire_dimensions (nc_data);

  /* Constrain to given to given dimensionality */
  t8_nc_constrain_dimensions_to_given_dimensionality (nc_data, dimensionality);

  /* Inquire the coordinate dimension data */
  t8_nc_inquire_dimensions_data (nc_data);

  /* Deallocate the netCDF specific variable data (which was allocated in 't8_nc_inquire_variables') and is not used anymore */
  t8_nc_deallocate_var_specific_data (nc_data);

  /* The mesh has to be built hereafter */
  /* .................................. */

#endif
}

void
t8_nc_construct_mesh_for_variables (t8_nc_data_t nc_data, const int num_variables, const char** const variable_names,
                                    const size_t* start_ptr, const size_t* count_ptr,
                                    const enum t8_nc_geo_mesh_type mesh_type, const enum t8_nc_geo_mesh_form mesh_form,
                                    const enum t8_nc_geo_mesh_elements mesh_elems)
{
#ifdef T8_WITH_NETCDF

  /* Inquire the coordinate-dimensions information */
  t8_nc_inquire_dimensions (nc_data);

  /* Inquire the coordinate-dimension's data */
  t8_nc_inquire_dimensions_data (nc_data);

  /* Inquire the variable's information */
  t8_nc_inquire_variables (nc_data, num_variables, variable_names);

  /* Constrain the dimensions based on the start and count arrays.
     * This needs to be done after the meta information about the variables have been gathered,
     * because it depends on the axis ordering of the variables */
  t8_nc_constrain_dimensions_to_given_hyperslab (nc_data, start_ptr, count_ptr);

  /* Inquire the variable's data */
  t8_nc_inquire_variables_data (nc_data, start_ptr, count_ptr);

  /* Deallocate the netCDF specific variable data (which was allocated in 't8_nc_inquire_variables') and is not used anymore */
  t8_nc_deallocate_var_specific_data (nc_data);

  /* Build the initial mesh */
  t8_nc_build_initial_mesh (nc_data, mesh_type, mesh_form, mesh_elems);

#endif
}

#ifdef T8_WITH_NETCDF
void
t8_netcdf_exit (const int _err_code, const char* _location)
{
  std::cout << "T8_NETCDF_EXIT is invoked..." << std::endl
            << "A netCDF-Error occurred, Code: " << _err_code << std::endl
            << nc_strerror (_err_code) << std::endl
            << "In: " << _location << std::endl;
  std::exit (EXIT_FAILURE);
}
#endif
