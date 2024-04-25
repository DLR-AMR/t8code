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
#include <t8_netcdf/t8_nc.hxx>
#include <t8_netcdf/t8_nc_geo_utilities.hxx>

#include <utility>
#include <array>
#include <vector>

void
t8_nc_exit (const int _err_code, const char* _location)
{
#ifdef T8_WITH_NETCDF
  t8_errorf ("A netCDF error occurred: Code %d at %s\n%s", _err_code, _location, nc_strerror (_err_code));
#endif
}

static int
nc_open_serial (const char* path_to_file)
{
#ifdef T8_WITH_NETCDF
  int err;
  int ncid;

  /* Open the file without explicit parallel access */
  err = nc__open (path_to_file, NC_NOWRITE, NULL, &ncid);
  t8_nc_check_error (err);

  return ncid;
#else
  t8_errorf ("t8code is not linked against netCDF. Please recompile the build with netCDF\n.");
  return T8_NC_ERR;
#endif
}

static int
nc_open_parallel (const char* path_to_file, sc_MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF_PAR
  int err;
  int ncid;

  int comm_size;
  err = sc_MPI_Comm_size (comm, &comm_size);
  SC_CHECK_MPI (err);

  MPI_Info info = MPI_INFO_NULL;
  err = nc_open_par (path_to_file, NC_NOWRITE, comm, info, &ncid);
  t8_nc_check_error (err);

  return ncid;
#else
  t8_errorf ("t8code is not linked against netCDF (at least no parallel access functions are available). Please "
             "recompile the build with netCDF\n.");
  return T8_NC_ERR;
#endif
}

void
t8_nc_data_t::nc_open (const std::string& path_to_file, const t8_nc_opening_mode_t mode, const MPI_Comm comm)
{
#ifdef T8_WITH_NETCDF
#ifdef T8_WITH_NETCDF_PAR
  if (mode == t8_nc_opening_mode_t::PARALLEL) {
    ncid_ = nc_open_parallel (path_to_file.c_str (), comm);
  }
  else

#endif /* T8_WITH_NETCDF_PAR */
  {
    if (mode == t8_nc_opening_mode_t::PARALLEL) {
      t8_debugf ("The netCDF file is ought to be opened for parallel access, although the parallel functionality is "
                 "not accessible.");
      t8_debugf ("The file is opened for serial access.");
    }
    ncid_ = nc_open_serial (path_to_file.c_str ());
  }

#else
  t8_errorf ("t8code is not linked against netCDF. Please recompile the build with netCDF\n.");
  return T8_NC_ERR;
#endif
}

void
t8_nc_data_t::set_hint_longitude_dimension (const int longitude_dimension_id)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (longitude_dimension_id >= 0);
  coordinate_dimension_ids_[t8_nc_dimension_t::LON] = longitude_dimension_id;
#endif
}

void
t8_nc_data_t::set_hint_latitude_dimension (const int latitude_dimension_id)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (latitude_dimension_id >= 0);
  coordinate_dimension_ids_[t8_nc_dimension_t::LAT] = latitude_dimension_id;
#endif
}

void
t8_nc_data_t::set_hint_height_dimension (const int height_dimension_id)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (height_dimension_id >= 0);
  coordinate_dimension_ids_[t8_nc_dimension_t::LEV] = height_dimension_id;
#endif
}

void
t8_nc_data_t::set_hint_time_dimension (const int time_dimension_id)
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (time_dimension_id >= 0);
  coordinate_dimension_ids_[t8_nc_dimension_t::TIME] = time_dimension_id;
#endif
}

/** Inquire the dimension lengths and ids of the (geo-spatial) coordinate dimensions as well as their variables ids */
void
t8_nc_data_t::inquire_coordinate_dimensions ()
{
#ifdef T8_WITH_NETCDF
  std::basic_string<char> nc_name;
  nc_name.reserve (NC_MAX_NAME + 1);

  size_t dim_length;

  /* Reserve memory for all 'dimension_lengths_' and 'dimension_names_' */
  dimension_lengths_.reserve (num_dimensions_);
  dimension_names_.reserve (num_dimensions_);

  const int& number_dimensions = num_dimensions_;

  /* Loop over all dimension within the netCDF file */
  for (int dim_index = 0; dim_index < number_dimensions; ++dim_index) {
    /* Get the name and the length of the dimension which corresponds to the current id */
    int err = nc_inq_dim (ncid_, dim_index, &nc_name[0], &dim_length);
    t8_nc_check_error (err);

    /* Check if a hint to the dimension id was given (concerning the coordinate dimensions) */
    if (coordinate_dimension_ids_[t8_nc_dimension_t::LON] == dim_index) {
      coordinate_lengths_[t8_nc_dimension_t::LON] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::LAT] == dim_index) {
      coordinate_lengths_[t8_nc_dimension_t::LAT] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::LEV] == dim_index) {
      coordinate_lengths_[t8_nc_dimension_t::LEV] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::TIME] == dim_index) {
      coordinate_lengths_[t8_nc_dimension_t::TIME] = static_cast<t8_gloidx_t> (dim_length);
    }
    /* In case no hints concerning the coordinate dimensions were given, we test for some default names */
    /* Check if the name is equal to any of the considered (geo-spatial) coordinate variables */
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::LON] == T8_NC_DIMENSION_NOT_CONSIDERED
             && (strcmp (nc_name.c_str (), "lon") == 0 || strcmp (nc_name.c_str (), "longitude") == 0)) {
      coordinate_dimension_ids_[t8_nc_dimension_t::LON] = dim_index;
      coordinate_lengths_[t8_nc_dimension_t::LON] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::LAT] == T8_NC_DIMENSION_NOT_CONSIDERED
             && (strcmp (nc_name.c_str (), "lat") == 0 || strcmp (nc_name.c_str (), "latitude") == 0)) {
      coordinate_dimension_ids_[t8_nc_dimension_t::LAT] = dim_index;
      coordinate_lengths_[t8_nc_dimension_t::LAT] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::LEV] == T8_NC_DIMENSION_NOT_CONSIDERED
             && (strcmp (nc_name.c_str (), "lev") == 0 || strcmp (nc_name.c_str (), "level") == 0
                 || strcmp (nc_name.c_str (), "height") == 0 || strcmp (nc_name.c_str (), "elevation") == 0)) {
      coordinate_dimension_ids_[t8_nc_dimension_t::LEV] = dim_index;
      coordinate_lengths_[t8_nc_dimension_t::LEV] = static_cast<t8_gloidx_t> (dim_length);
    }
    else if (coordinate_dimension_ids_[t8_nc_dimension_t::TIME] == T8_NC_DIMENSION_NOT_CONSIDERED
             && strcmp (nc_name.c_str (), "time") == 0) {
      coordinate_dimension_ids_[t8_nc_dimension_t::TIME] = dim_index;
      coordinate_lengths_[t8_nc_dimension_t::TIME] = static_cast<t8_gloidx_t> (dim_length);
    }

    /* Save the name and the length of each dimension within the netCDF file */
    dimension_names_.emplace_back (nc_name.c_str ());
    dimension_lengths_.push_back (static_cast<t8_gloidx_t> (dim_length));
  }

  //cmc_debug_msg("The following geo-spatial coordinate dimensions were retrieved from the netCDF file:");
  //cmc_debug_msg("LON: ", coordinate_lengths_[t8_nc_dimension_t::LON], ", Lat: ", coordinate_lengths_[t8_nc_dimension_t::LAT], ", LEV: ", coordinate_lengths_[t8_nc_dimension_t::LEV], ", T: ", coordinate_lengths_[t8_nc_dimension_t::TIME]);

#endif
}

void
t8_nc_data_t::inquire_coordinates ()
{
#ifdef T8_WITH_NETCDF
  /* Inquire the most general information about the supplied netCDF file (like number of diemnsion, number of global attributes, number of variables, etc.)*/
  int err = nc_inq (ncid_, &num_dimensions_, NULL, &num_global_attributes_, &id_unlimited_dimension_);
  t8_nc_check_error (err);

  //TODO: When the coordinates are read in. They need to be adjusted to the given hyperslab later
  inquire_coordinate_dimensions ();
  //InquireCoordinateVariables();
  //InquireCoordinateData();

#endif
}

static int
get_variable_id (const int ncid, const std::string& variable_name)
{
#ifdef T8_WITH_NETCDF
  int variable_id;

  /* Inquire the netCDF internal variable id by the given variable's name */
  const int err = nc_inq_varid (ncid, variable_name.c_str (), &variable_id);
  t8_nc_check_error (err);

  return variable_id;
#else
  return T8_NC_ERR;
#endif
}

static std::pair<bool, t8_universal_type_t>
check_variable_for_attribute (const int ncid, const int variable_id, const std::string& attribute_name)
{
#ifdef T8_WITH_NETCDF
  int attribute_type;

  int err = nc_inq_atttype (ncid, variable_id, attribute_name.c_str (), &attribute_type);

  t8_universal_type_t attribute_value;

  if (err != NC_NOERR) {
    /* The specific attribute has not been found */
    return std::make_pair (false, attribute_value);
  }

  //cmc_debug_msg("The attribute ", attribute_name, " has been found.");

  /* Otherwise, if the element has been found, we check it's type and inquire the attribute */
  switch (attribute_type) {
  case NC_DOUBLE:
    double double_att;
    err = nc_get_att_double (ncid, variable_id, attribute_name.c_str (), &double_att);
    t8_nc_check_error (err);
    attribute_value = double_att;
    break;
  case NC_INT:
    int32_t int_att;
    err = nc_get_att_int (ncid, variable_id, attribute_name.c_str (), &int_att);
    t8_nc_check_error (err);
    attribute_value = int_att;
    break;
  case NC_FLOAT:
    float float_att;
    err = nc_get_att_float (ncid, variable_id, attribute_name.c_str (), &float_att);
    t8_nc_check_error (err);
    attribute_value = float_att;
    break;
  case NC_SHORT:
    int16_t short_att;
    err = nc_get_att_short (ncid, variable_id, attribute_name.c_str (), &short_att);
    t8_nc_check_error (err);
    attribute_value = short_att;
    break;
  case NC_USHORT:
    uint16_t ushort_att;
    err = nc_get_att_ushort (ncid, variable_id, attribute_name.c_str (), &ushort_att);
    t8_nc_check_error (err);
    attribute_value = ushort_att;
    break;
  case NC_UINT:
    uint32_t uint_att;
    err = nc_get_att_uint (ncid, variable_id, attribute_name.c_str (), &uint_att);
    t8_nc_check_error (err);
    attribute_value = uint_att;
    break;
  case NC_INT64:
    long long int longlong_att;
    err = nc_get_att_longlong (ncid, variable_id, attribute_name.c_str (), &longlong_att);
    t8_nc_check_error (err);
    attribute_value = static_cast<int64_t> (longlong_att);
    break;
  case NC_UINT64:
    unsigned long long int ulonglong_att;
    err = nc_get_att_ulonglong (ncid, variable_id, attribute_name.c_str (), &ulonglong_att);
    t8_nc_check_error (err);
    attribute_value = static_cast<uint64_t> (ulonglong_att);
    break;
  default:
    t8_debugf ("The given attribute could not be inquired, since its data type is not supported.");
    //cmc_debug_msg("The attribute ", attribute_name, " of variable with ID: ", variable_id, " could not be inquired (and will not be applied), since its data type is not supported.");
  }

  return std::make_pair (true, std::move (attribute_value));
#else
  return std::make_pair (false, T8_NC_ERR);
#endif
}

static t8_nc_data_layout_t
get_layout_from_dimension_ordering (const std::vector<t8_nc_dimension_t>& dimension_ordering)
{
  if (dimension_ordering.size () == 2) {
    /* Two dimensional layout */
    if (dimension_ordering[0] == t8_nc_dimension_t::LON && dimension_ordering[1] == t8_nc_dimension_t::LAT)
      return t8_nc_data_layout_t::LON_LAT;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LON && dimension_ordering[1] == t8_nc_dimension_t::LEV)
      return t8_nc_data_layout_t::LON_LEV;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LAT && dimension_ordering[1] == t8_nc_dimension_t::LON)
      return t8_nc_data_layout_t::LAT_LON;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LAT && dimension_ordering[1] == t8_nc_dimension_t::LEV)
      return t8_nc_data_layout_t::LAT_LEV;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LEV && dimension_ordering[1] == t8_nc_dimension_t::LON)
      return t8_nc_data_layout_t::LEV_LON;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LEV && dimension_ordering[1] == t8_nc_dimension_t::LAT)
      return t8_nc_data_layout_t::LEV_LAT;
    else {
      t8_errorf ("The two-dimensional data layout could not be deduced from the dimension ordering.");
      return t8_nc_data_layout_t::LAYOUT_UNDEFINED;
    }
  }
  else if (dimension_ordering.size () == 3) {
    /* Three dimensional layout */
    if (dimension_ordering[0] == t8_nc_dimension_t::LEV && dimension_ordering[1] == t8_nc_dimension_t::LON
        && dimension_ordering[2] == t8_nc_dimension_t::LAT)
      return t8_nc_data_layout_t::LEV_LON_LAT;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LEV && dimension_ordering[1] == t8_nc_dimension_t::LAT
             && dimension_ordering[2] == t8_nc_dimension_t::LON)
      return t8_nc_data_layout_t::LEV_LAT_LON;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LAT && dimension_ordering[1] == t8_nc_dimension_t::LEV
             && dimension_ordering[2] == t8_nc_dimension_t::LON)
      return t8_nc_data_layout_t::LAT_LEV_LON;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LAT && dimension_ordering[1] == t8_nc_dimension_t::LON
             && dimension_ordering[2] == t8_nc_dimension_t::LEV)
      return t8_nc_data_layout_t::LAT_LON_LEV;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LON && dimension_ordering[1] == t8_nc_dimension_t::LAT
             && dimension_ordering[2] == t8_nc_dimension_t::LEV)
      return t8_nc_data_layout_t::LON_LAT_LEV;
    else if (dimension_ordering[0] == t8_nc_dimension_t::LON && dimension_ordering[1] == t8_nc_dimension_t::LEV
             && dimension_ordering[2] == t8_nc_dimension_t::LAT)
      return t8_nc_data_layout_t::LON_LEV_LAT;
    else {
      t8_errorf ("The three-dimensional data layout could not be deduced from the dimension ordering.");
      return t8_nc_data_layout_t::LAYOUT_UNDEFINED;
    }
  }
  else {
    t8_errorf ("Only 2D and 3D variables are currently supported, therefore the t8_nc_data_layout_t could not be "
               "deduced from the dimension ordering.");
    return t8_nc_data_layout_t::LAYOUT_UNDEFINED;
  }
}

static t8_nc_data_layout_t
get_data_layout (const int dimensionality, const t8_nc_hyperslab_t& global_hyperslab,
                 const t8_nc_coordinate_array_t<int>& coordinate_dimension_ids_,
                 const std::array<int, NC_MAX_VAR_DIMS>& variable_dimension_ids)
{
  std::vector<t8_nc_dimension_t> data_dimensions;
  data_dimensions.reserve (dimensionality);

  int found_dimensions = 0;

  for (auto dim_iter = variable_dimension_ids.begin ();
       found_dimensions < dimensionality && dim_iter != variable_dimension_ids.end (); ++dim_iter) {
    if (const auto coordinate_dim_id_iter
        = std::find_if (coordinate_dimension_ids_.begin (), coordinate_dimension_ids_.end (),
                        [&] (const int& arg) { return arg == *dim_iter; });
        coordinate_dim_id_iter != coordinate_dimension_ids_.end ()) {
      /* If the current dimension id corresponds to a coordinate dimension id (which is considered), we add the dimension to the vector */
      const t8_nc_dimension_t dimension
        = static_cast<t8_nc_dimension_t> (std::distance (coordinate_dimension_ids_.begin (), coordinate_dim_id_iter));
      if (global_hyperslab.get_dimension_length (dimension) > 1) {
        data_dimensions.push_back (dimension);
        ++found_dimensions;
      }
    }
  }

  if (found_dimensions < dimensionality) {
    t8_errorf ("The t8_nc_data_layout_t could not be deduced. The variable is not only dependent on geo-spatial "
               "coordinate dimensions");
  }

  return get_layout_from_dimension_ordering (data_dimensions);
}

#if 0


static
std::vector<t8_nc_hyperslab_t>
determine_data_distribution(const t8_nc_hyperslab_t& global_hyperslab)
{
    //TODO:: Update for the parallel case
    return std::vector<t8_nc_hyperslab_t>{global_hyperslab};
}

static
std::pair<std::vector<size_t>, std::vector<size_t>>
get_start_and_count_values_for_variable(const t8_nc_hyperslab_t& hyperslab, const t8_nc_coordinate_array_t<int>& coordinate_dimension_ids, int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    std::vector<size_t> start_vals(num_dimensions, 0);
    std::vector<size_t> count_vals(num_dimensions, 1);

    for (int iter = 0; iter < num_dimensions; ++iter)
    {
        if (const auto coordinate_dim_id_iter = std::find_if(coordinate_dimension_ids.begin(), coordinate_dimension_ids.end(), [&](const int& arg){return arg == dimension_ids[iter];});
            coordinate_dim_id_iter != coordinate_dimension_ids.end())
        {
            const t8_nc_dimension_t dimension = static_cast<t8_nc_dimension_t>(std::distance(coordinate_dimension_ids.begin(), coordinate_dim_id_iter));
            start_vals[iter] = static_cast<size_t>(hyperslab.get_dimension_start(dimension));
            count_vals[iter] = static_cast<size_t>(hyperslab.get_dimension_length(dimension));
        }
        
    }

    return std::make_pair(std::move(start_vals), std::move(count_vals));
}

template<typename T>
static InputVar
SetUpInputVariable(const int ncid, const t8_nc_coordinate_array_t<int>& coordinate_dimension_ids, const int variable_id, std::string&& variable_name, const t8_nc_data_layout_t layout,
                   const t8_gloidx_t num_values, std::vector<t8_nc_hyperslab_t>&& hyperslabs, GeoDomain&& global_domain,
                   const int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    InputVariable<T> variable(std::move(variable_name), variable_id, layout);

    variable.SetGlobalDomain(std::move(global_domain));

    std::vector<T> local_data(num_values);
    t8_gloidx_t offset = 0;
    for (auto hs_iter = hyperslabs.begin(); hs_iter != hyperslabs.end(); ++hs_iter)
    {
        const auto [start_vals, count_vals] = get_start_and_count_values_for_variable(*hs_iter, coordinate_dimension_ids, num_dimensions, dimension_ids);

        const int err = nc_get_vara(ncid, variable_id, start_vals.data(), count_vals.data(), static_cast<void*>(&local_data[offset]));
        t8_nc_check_error(err);

        offset += hs_iter->get_number_coordinates();
    }
    
    variable.SetDataAndCoordinates(std::move(local_data), std::move(hyperslabs));

    return InputVar(std::move(variable));
}

InputVar
t8_nc_data_t::SetupVariableData(const int variable_nc_type, const int variable_id, std::string&& variable_name,
                          const t8_nc_data_layout_t layout, const t8_gloidx_t num_values, std::vector<t8_nc_hyperslab_t>&& hyperslabs,
                          GeoDomain&& global_domain, const int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    switch (variable_nc_type)
    {
        case NC_DOUBLE:
            return SetUpInputVariable<double>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_FLOAT:
            return SetUpInputVariable<float>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_INT:
            return SetUpInputVariable<int32_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UINT:
            return SetUpInputVariable<uint32_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_INT64:
            return SetUpInputVariable<int64_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UINT64:
            return SetUpInputVariable<uint64_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_SHORT:
            return SetUpInputVariable<int16_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_USHORT:
            return SetUpInputVariable<uint16_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_BYTE:
            return SetUpInputVariable<int8_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UBYTE:
            return SetUpInputVariable<uint8_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_CHAR:
            return SetUpInputVariable<char>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        default:
            t8_errorf("The variable of type (NC-Type: ", variable_nc_type, ") is not supported.");
            return InputVar();
    }
}


InputVar
t8_nc_data_t::InquireVariable(const t8_nc_hyperslab_t& hyperslab, std::string&& variable_name)
{
    const int variable_id = GetVariableId(ncid_, variable_name);

    int variable_nc_type;
    int num_dimensions;
    std::array<int, NC_MAX_VAR_DIMS> dimension_ids;

    /* Inquire some information about the variable, e.g. the data type and its dimensions */
    int err = nc_inq_var(ncid_, variable_id, NULL, &variable_nc_type, &num_dimensions, dimension_ids.data(), NULL);
    t8_nc_check_error(err);

    /* Check some attributes of the variable */
    const auto [is_missing_value_present, missing_value] = check_variable_for_attribute(ncid_, variable_id, "missing_value");
    const auto [is_add_offset_present, add_offset] = check_variable_for_attribute(ncid_, variable_id, "add_offset");
    const auto [is_scale_factor_present, scale_factor] = check_variable_for_attribute(ncid_, variable_id, "scale_factor");

    const int data_dimensionality = hyperslab.GetDimensionality();

    /* Determine the data layout of the variable */
    const t8_nc_data_layout_t data_layout = get_data_layout(data_dimensionality, hyperslab, coordinate_dimension_ids_, dimension_ids);

    /* TODO: Add parallelization */
    //for now in serial only

    std::vector<t8_nc_hyperslab_t> local_hyperslabs = determine_data_distribution(hyperslab);

    /* Get the number of values we will read from the file (for memory allocation) */
    t8_gloidx_t num_local_data_values = 0;
    for (auto hs_iter = local_hyperslabs.begin(); hs_iter != local_hyperslabs.end(); ++hs_iter)
    {
        num_local_data_values += hs_iter->GetNumberCoordinates();
    }

    /* Inquire the data of the variable and construct a InputVariable with it */
    InputVar variable = SetupVariableData(variable_nc_type, variable_id, std::move(variable_name), data_layout, num_local_data_values, std::move(local_hyperslabs), TransformHyperslabToGeoDomain(hyperslab), num_dimensions, dimension_ids);

    if (is_missing_value_present)
        variable.SetMissingValue(missing_value);
    
    if (is_add_offset_present)
        variable.SetAddOffset(add_offset);
    
    if (is_scale_factor_present)
        variable.SetScaleFactor(scale_factor);

    return variable;
}

#endif

//TODO: Needs to be updated for the real input variable
InputVar
t8_nc_data_t::inquire_variable (const t8_nc_hyperslab_t& hyperslab, std::string&& variable_name)
{
#ifdef T8_WITH_NETCDF
  const int variable_id = get_variable_id (ncid_, variable_name);

  int variable_nc_type;
  int num_dimensions;
  std::array<int, NC_MAX_VAR_DIMS> dimension_ids;

  /* Inquire some information about the variable, e.g. the data type and its dimensions */
  int err = nc_inq_var (ncid_, variable_id, NULL, &variable_nc_type, &num_dimensions, dimension_ids.data (), NULL);
  t8_nc_check_error (err);

  /* Check some attributes of the variable */
  const auto [is_missing_value_present, missing_value]
    = check_variable_for_attribute (ncid_, variable_id, "missing_value");
  const auto [is_add_offset_present, add_offset] = check_variable_for_attribute (ncid_, variable_id, "add_offset");
  const auto [is_scale_factor_present, scale_factor]
    = check_variable_for_attribute (ncid_, variable_id, "scale_factor");

  const int data_dimensionality = hyperslab.get_dimensionality ();

  /* Determine the data layout of the variable */
  const t8_nc_data_layout_t data_layout
    = get_data_layout (data_dimensionality, hyperslab, coordinate_dimension_ids_, dimension_ids);

  InputVar variable;
  variable.variable_id_ = variable_id;
  variable.dimensionality_ = data_dimensionality;
  variable.initial_layout_ = data_layout;
  variable.global_domain_ = transform_hyperslab_to_geo_domain (hyperslab);

  if (is_missing_value_present)
    variable.missing_value_ = missing_value;

  if (is_add_offset_present)
    variable.add_offset_ = add_offset;

  if (is_scale_factor_present)
    variable.scale_factor_ = scale_factor;

  return variable;
#else
  return InputVar ();
#endif
}

[[nodiscard]] std::vector<InputVar>&&
t8_nc_data_t::transfer_data ()
{
  if (_data_has_been_transferred_)
    t8_errorf ("The data already has been transferred.");

  return std::move (variables_);
}
