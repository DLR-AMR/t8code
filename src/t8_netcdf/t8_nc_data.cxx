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

#include <t8_netcdf/t8_nc_data.hxx>
#include <vector>
#include <t8.h>

#if 0
/**
 * \brief A struct describing and holding geo-spatial data
 * 
 */
struct t8_geo_var
{
 private:
 public:
  t8_geo_var () {};
  t8_geo_var (std::string _name): name { _name } {};
  t8_geo_var (std::string _name, const t8_geo_data_type type): name { _name }, type { type } {};
  t8_geo_var (const t8_geo_data_type type, const size_t num_elements): type { type }
  {
    data = sc_array_new_count (t8_nc_type_to_bytes (type), num_elements);
  };
  t8_geo_var (std::string _name, const t8_geo_data_type type, const size_t num_elements): name { _name }, type { type }
  {
    data = sc_array_new_count (t8_nc_type_to_bytes (type), num_elements);
  };
  ~t8_geo_var ()
  {
    if (data != nullptr) {
      sc_array_destroy (data);
    }
    if (data_new != nullptr) {
      sc_array_destroy (data_new);
    }
  };

  /********** General Variable Data **********/
  /* Meta information about the variable */
  std::string name {};  //!< The name of the variable
  t8_geo_data_type type {
    t8_geo_data_type::T8_GEO_DATA_UNDEFINED
  };  //!< The data type of this variable (e.g. T8_INT32_T for int32_t data)

  /* A pointer to the data array of the variable */
  sc_array_t* data { nullptr };  //!< The 'actual' data of the variable
  sc_array_t* data_new {
    nullptr
  };  //!< A convenience variable which may be used if the data is transformed/adapted/partitioned etc.

  t8_nc_data_ordering data_ordering {
    t8_nc_data_ordering::T8_LAYOUT_UNDEFINED
  };  //!< An enumerator indicating the layout and ordering of the data

  t8_universal_type_t
    missing_value;  //!< A value indicating that no 'real' value exists at this position (this value is just a placeholder for a missing value)
  t8_universal_type_t add_offset { 0 };    //!< Indicating the offset of the data if it was shifted
  t8_universal_type_t scale_factor { 1 };  //!< Indicating the scaling of the data

  bool missing_values_present_within_the_data {
    false
  };  //!< Flag indicating whether or not missing values exists within the data
  bool applied_offset_and_scaling {
    false
  };  //!< Flag indicating whether the offset and scaling has been applied to the data or not

  void* user_data { nullptr };  //!< A member variable capable of holding user_data
};

#ifdef T8_WITH_NETCDF
t8_geo_data_type
t8_nc_geo_variable_nc_type_to_t8_geo_data_type (nc_type nc_data_type)
{
  /** The custom 't8_geo_data_type' has the same layout/specifications as the nc_type from netCDF 
     *  Therefore, it suffices to just cast the value */
  return static_cast<t8_geo_data_type> (nc_data_type);
}
#endif

size_t
t8_nc_type_to_bytes (const enum t8_geo_data_type type)
{
  T8_ASSERT (type >= 0 && type < t8_geo_data_type::T8_GEO_DATA_NUM_TYPES);
  constexpr std::array<size_t, t8_geo_data_type::T8_GEO_DATA_NUM_TYPES> bytes {
#ifdef _cpp_lib_byte
    sizeof (std::byte)
#else
    sizeof (unsigned char)
#endif
      ,
    sizeof (int8_t),
    sizeof (char),
    sizeof (int16_t),
    sizeof (int32_t),
    sizeof (float),
    sizeof (double),
    sizeof (uint8_t),
    sizeof (uint16_t),
    sizeof (uint32_t),
    sizeof (int64_t),
    sizeof (uint64_t)
  };

  return bytes[type];
}

t8_geo_var_t
t8_nc_create_geo_variable (const char* name)
{
  /* Allocate a new variable and return a pointer to it */
  return new t8_geo_var (std::string (name));
}

t8_geo_var_t
t8_nc_create_geo_variable (const char* name, const enum t8_geo_data_type type)
{
  /* Allocate a new variable and return a pointer to it */
  return new t8_geo_var (std::string (name), type);
}

t8_geo_var_t
t8_nc_create_geo_variable (const char* name, const enum t8_geo_data_type type, const size_t num_elements)
{
  /* Allocate a new variable and return a pointer to it */
  return new t8_geo_var (std::string (name), type, num_elements);
}

void
t8_nc_destroy_geo_variable (t8_geo_var_t var)
{
  /* Deallocate the variable */
  delete var;
}

void
t8_nc_geo_variable_allocate_data (t8_geo_var_t var, const t8_geo_data_type type, const size_t num_elements)
{
  /* Deallocate the previous array if there is one present */
  if (var->data != nullptr) {
    sc_array_destroy (var->data);
  }
  /* Allocate a new sc_array */
  var->data = sc_array_new_count (t8_nc_type_to_bytes (type), num_elements);
  /* Store the supplied data type */
  var->type = type;
}

t8_geo_data_type
t8_nc_geo_variable_get_type (t8_geo_var_t var)
{
  /* Return the type of the variable */
  return var->type;
}

size_t
t8_nc_geo_variable_get_num_elements (t8_geo_var_t var)
{
  /* Return the number of elements of this variable */
  return var->data->elem_count;
}

void*
t8_nc_geo_variable_get_data_ptr (t8_geo_var_t var)
{
  /* Return the underlying data pointer, if there is one */
  return static_cast<void*> (var->data != nullptr ? var->data->array : nullptr);
}

void
t8_nc_geo_variable_set_user_data (t8_geo_var_t var, void* data)
{
  /* Set the user data */
  var->user_data = data;
}

void*
t8_nc_geo_variable_get_user_data (t8_geo_var_t var)
{
  /* Get the user data */
  return var->user_data;
}

const char*
t8_nc_geo_variable_get_name (t8_geo_var_t var)
{
  /* Return the name of the variable */
  return var->name.c_str ();
}

void
t8_nc_geo_variable_set_name (t8_geo_var_t var, const char* name)
{
  /* Set the name of the variable */
  var->name = std::string (name);
}

t8_universal_type_t
t8_nc_geo_variable_get_missing_value (t8_geo_var_t var)
{
  return var->missing_value;
}

void
t8_nc_geo_variable_set_missing_value (t8_geo_var_t var, t8_universal_type_t missing_value)
{
  var->missing_value = missing_value;
}

t8_universal_type_t
t8_nc_geo_variable_get_scale_factor (t8_geo_var_t var)
{
  return var->scale_factor;
}

void
t8_nc_geo_variable_set_scale_factor (t8_geo_var_t var, t8_universal_type_t scale_factor)
{
  var->scale_factor = scale_factor;
}

t8_universal_type_t
t8_nc_geo_variable_get_add_offset (t8_geo_var_t var)
{
  return var->add_offset;
}

void
t8_nc_geo_variable_set_add_offset (t8_geo_var_t var, t8_universal_type_t offset)
{
  var->add_offset = offset;
}

static void
t8_nc_geo_variable_switch_data (t8_geo_var_t var)
{
  if (var->data != nullptr) {
    /* Deallocate the array */
    sc_array_destroy (var->data);
  }

  /* Exchange the pointers to the sc_arrays */
  var->data = var->data_new;

  /* Nullify the helper array pointer */
  var->data_new = nullptr;
}

void
t8_nc_geo_variable_crop_data_to_selection (t8_geo_var_t var, const size_t start_index, const size_t end_index)
{
  T8_ASSERT (var->data->elem_count > end_index && start_index <= end_index);

  /* Get the type of the variable */
  const enum t8_geo_data_type type = t8_nc_geo_variable_get_type (var);

  /* Allocate the data_new array for the cropped data */
  var->data_new = sc_array_new_count (t8_nc_type_to_bytes (type), end_index + 1 - start_index);

  /* Copy the data from the previous array to the cropped array */
  memcpy (var->data_new->array, var->data->array + t8_nc_type_to_bytes (type) * start_index,
          (end_index + 1 - start_index) * t8_nc_type_to_bytes (type));

  /* Switch both variable internal array, such that only the cropped one remains */
  t8_nc_geo_variable_switch_data (var);
}

void
t8_nc_geo_variable_set_data_ordering_scheme (t8_geo_var_t var, const t8_nc_data_ordering data_ordering)
{
  var->data_ordering = data_ordering;
}

t8_nc_data_ordering
t8_nc_geo_variable_get_data_ordering_scheme (t8_geo_var_t var)
{
  return var->data_ordering;
}
#endif