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

#ifndef T8_NC_DATA_HXX
#define T8_NC_DATA_HXX

#include <t8.h>

#ifdef T8_WITH_NETCDF
#include <netcdf.h>
#endif

#if 0
/* Enum for variable data types */
enum t8_geo_data_type {
  T8_GEO_DATA_UNDEFINED = -1,
  T8_BYTE,
  T8_INT8_T,
  T8_CHAR,
  T8_INT16_T,
  T8_INT32_T,
  T8_FLOAT,
  T8_DOUBLE,
  T8_UINT8_T,
  T8_UINT16_T,
  T8_UINT32_T,
  T8_INT64_T,
  T8_UINT64_T,
  T8_GEO_DATA_NUM_TYPES
};
#endif

#if 0
/* Enum for the data layout of the variable */
enum t8_geo_data_layout {
  T8_GEO_DATA_LAYOUT_UNDEFINED = -1,
  T8_GEO_DATA_SFC_ORDERED,
  T8_GEO_DATA_BLOCKED,
  T8_GEO_DATA_POINTS
};

/* An enumerator for data ordering */
enum t8_nc_data_ordering {
  T8_LAYOUT_UNDEFINED = -1,
  T8_2D_LAT_LON,
  T8_2D_LON_LAT,
  T8_2D_LAT_LEV,
  T8_2D_LEV_LAT,
  T8_2D_LON_LEV,
  T8_2D_LEV_LON,
  _INTERN_ID_END_CARTESIAN_2D,
  T8_3D_LAT_LON_LEV,
  T8_3D_LAT_LEV_LON,
  T8_3D_LEV_LAT_LON,
  T8_3D_LEV_LON_LAT,
  T8_3D_LON_LEV_LAT,
  T8_3D_LON_LAT_LEV,
  _INTERN_ID_END_CARTESIAN_3D,
  T8_2D_SFC,
  T8_3D_SFC
};

/* Pointer typedef for geo-spatial variables */
typedef struct t8_geo_var* t8_geo_var_t;

/* Define a type capable of holding different data types */
#ifdef __cpp_lib_variant
/* A typedef for a variant capable of holding different data types (if a variant is accessible) */
#include <variant>
typedef std::variant<unsigned char, int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t,
                     uint64_t>
  t8_universal_type_t;
#else
/* A variant equivalent formulation for compliance with C++ standards before C++17 */
typedef union t8_universal_type {
  unsigned char b;
  int8_t i8;
  char c;
  int16_t i16;
  int32_t i32;
  float f;
  double d;
  uint8_t ui8;
  uint16_t ui16;
  uint32_t ui32;
  int64_t i64;
  uint64_t ui64;
#ifdef __cplusplus
  /* Overload functions for convenience */
  template <typename T>
  auto
  operator= (T&& value) -> typename std::enable_if<std::is_arithmetic<T>::value, t8_universal_type&>::type;
#endif
} t8_universal_type_t;
#endif

size_t
t8_nc_type_to_bytes (const enum t8_geo_data_type type);

t8_geo_var_t
t8_nc_create_geo_variable (const char* name);

t8_geo_var_t
t8_nc_create_geo_variable (const char* name, const enum t8_geo_data_type type);

t8_geo_var_t
t8_nc_create_geo_variable (const char* name, const enum t8_geo_data_type type, const size_t num_elements);

void
t8_nc_destroy_geo_variable (t8_geo_var_t var);

void*
t8_nc_geo_variable_get_data_ptr (t8_geo_var_t var);

void
t8_nc_geo_variable_set_user_data (t8_geo_var_t var, void* data);

void*
t8_nc_geo_variable_get_user_data (t8_geo_var_t var);

const char*
t8_nc_geo_variable_get_name (t8_geo_var_t var);

void
t8_nc_geo_variable_set_name (t8_geo_var_t var, const char* name);

t8_universal_type_t
t8_nc_geo_variable_get_missing_value (t8_geo_var_t var);

void
t8_nc_geo_variable_set_missing_value (t8_geo_var_t var, t8_universal_type_t missing_value);

t8_universal_type_t
t8_nc_geo_variable_get_scale_factor (t8_geo_var_t var);

void
t8_nc_geo_variable_set_scale_factor (t8_geo_var_t var, t8_universal_type_t scale_factor);

t8_universal_type_t
t8_nc_geo_variable_get_add_offset (t8_geo_var_t var);

void
t8_nc_geo_variable_set_add_offset (t8_geo_var_t var, t8_universal_type_t offset);

t8_geo_data_type
t8_nc_geo_variable_get_type (t8_geo_var_t var);

size_t
t8_nc_geo_variable_get_num_elements (t8_geo_var_t var);

void
t8_nc_geo_variable_crop_data_to_selection (t8_geo_var_t var, const size_t start_index, const size_t end_index);

void
t8_nc_geo_variable_allocate_data (t8_geo_var_t var, const t8_geo_data_type type, const size_t num_elements);

void
t8_nc_geo_variable_set_data_ordering_scheme (t8_geo_var_t var, const t8_nc_data_ordering data_ordering);

t8_nc_data_ordering
t8_nc_geo_variable_get_data_ordering_scheme (t8_geo_var_t var);

#ifdef T8_WITH_NETCDF
t8_geo_data_type
t8_nc_geo_variable_nc_type_to_t8_geo_data_type (nc_type nc_data_type);
#endif

#if defined __cplusplus && !defined __cpp_lib_variant
template <typename T, typename U>
constexpr bool t8_is_decay_same = std::is_same<std::decay_t<T>, U>::value;

template <typename T>
t8_geo_data_type
t8_decay_to_geo_data_type (T)
{
  if (t8_is_decay_same<T, unsigned char>) {
    return t8_geo_data_type::T8_BYTE;
  }
  else if (t8_is_decay_same<T, int8_t>) {
    return t8_geo_data_type::T8_INT8_T;
  }
  else if (t8_is_decay_same<T, char>) {
    return t8_geo_data_type::T8_CHAR;
  }
  else if (t8_is_decay_same<T, int16_t>) {
    return t8_geo_data_type::T8_INT16_T;
  }
  else if (t8_is_decay_same<T, int32_t>) {
    return t8_geo_data_type::T8_INT32_T;
  }
  else if (t8_is_decay_same<T, float>) {
    return t8_geo_data_type::T8_FLOAT;
  }
  else if (t8_is_decay_same<T, double>) {
    return t8_geo_data_type::T8_DOUBLE;
  }
  else if (t8_is_decay_same<T, uint8_t>) {
    return t8_geo_data_type::T8_UINT8_T;
  }
  else if (t8_is_decay_same<T, uint16_t>) {
    return t8_geo_data_type::T8_UINT16_T;
  }
  else if (t8_is_decay_same<T, uint32_t>) {
    return t8_geo_data_type::T8_UINT32_T;
  }
  else if (t8_is_decay_same<T, int64_t>) {
    return t8_geo_data_type::T8_INT64_T;
  }
  else if (t8_is_decay_same<T, uint64_t>) {
    return t8_geo_data_type::T8_UINT64_T;
  }
  else {
    t8_errorf ("The supplied data type does not match any of the types of t8_universal_type_t");
    return t8_geo_data_type::T8_GEO_DATA_UNDEFINED;
  }
}

template <typename T>
auto
t8_universal_type::operator= (T&& value) ->
  typename std::enable_if<std::is_arithmetic<T>::value, t8_universal_type&>::type
{
  t8_debugf ("Custom union assign\n");
  /* Declare an universal type */
  t8_universal_type_t& t = *this;

  /* Assign a proper value*/
  switch (t8_decay_to_geo_data_type (T ())) {
  case t8_geo_data_type::T8_BYTE:
    t.b = static_cast<unsigned char> (value);
    break;
  case t8_geo_data_type::T8_INT8_T:
    t.i8 = static_cast<int8_t> (value);
    break;
  case t8_geo_data_type::T8_CHAR:
    t.c = static_cast<char> (value);
    break;
  case t8_geo_data_type::T8_INT16_T:
    t.i16 = static_cast<int16_t> (value);
    break;
  case t8_geo_data_type::T8_INT32_T:
    t.i32 = static_cast<int8_t> (value);
    break;
  case t8_geo_data_type::T8_FLOAT:
    t.f = static_cast<float> (value);
    break;
  case t8_geo_data_type::T8_DOUBLE:
    t.d = static_cast<double> (value);
    break;
  case t8_geo_data_type::T8_UINT8_T:
    t.ui8 = static_cast<uint8_t> (value);
    break;
  case t8_geo_data_type::T8_UINT16_T:
    t.ui16 = static_cast<uint16_t> (value);
    break;
  case t8_geo_data_type::T8_UINT32_T:
    t.ui32 = static_cast<uint32_t> (value);
    break;
  case t8_geo_data_type::T8_INT64_T:
    t.i64 = static_cast<int64_t> (value);
    break;
  case t8_geo_data_type::T8_UINT64_T:
    t.ui64 = static_cast<uint64_t> (value);
    break;
  default:
    t8_errorf ("The supplied data type does not match any of the types of t8_universal_type_t");
  }

  return *this;
}
#endif

#endif

#endif /* !T8_NC_DATA_HXX */
