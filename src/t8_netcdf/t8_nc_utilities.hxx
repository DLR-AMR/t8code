/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

/**
 * \file t8_nc_utilities.hxx This file collects some utility functions for convenient use within the t8code-netCDF functionality
 */
#ifndef T8_NC_UTILITIES_HXX
#define T8_NC_UTILITIES_HXX

#include <t8.h>
#include <t8_netcdf/t8_nc_dimension_interval.hxx>

#include <type_traits>
#include <variant>

/* A typedef for a variant capable of holding different data types (if a variant is accessible) */
typedef std::variant<unsigned char, int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t,
                     uint64_t>
  t8_universal_type_t;

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

template <typename T>
class t8_nc_coordinate_array_t {
 public:
  t8_nc_coordinate_array_t () = default;
  t8_nc_coordinate_array_t (const T fill_value): coordinates { fill_value, fill_value, fill_value, fill_value } {};

  T&
  operator[] (int idx)
  {
    T8_ASSERT (idx < t8_nc_dimension_t::NUM_COORDINATES);
    return coordinates[idx];
  }
  const T&
  operator[] (int idx) const
  {
    T8_ASSERT (idx < t8_nc_dimension_t::NUM_COORDINATES);
    return coordinates[idx];
  }

  typename std::array<T, t8_nc_dimension_t::NUM_COORDINATES>::iterator
  begin ()
  {
    return coordinates.begin ();
  }
  typename std::array<T, t8_nc_dimension_t::NUM_COORDINATES>::iterator
  end ()
  {
    return coordinates.end ();
  }

  typename std::array<T, t8_nc_dimension_t::NUM_COORDINATES>::const_iterator
  begin () const
  {
    return coordinates.begin ();
  }
  typename std::array<T, t8_nc_dimension_t::NUM_COORDINATES>::const_iterator
  end () const
  {
    return coordinates.end ();
  }

  std::array<T, t8_nc_dimension_t::NUM_COORDINATES> coordinates { 0, 0, 0, 0 };
};

/**
 * \brief Using integer exponentiation by squaring
 * 
 * \param base The base value to be exponentiated
 * \param exponent The exponent
 * \return int The resulting value is the integer resulting from (base)^(exponent)
 * TODO: maybe check for overflow during computation
 */
template <typename Integer>
auto
int_pow (Integer base, Integer exponent) -> typename std::enable_if<std::is_integral<Integer>::value, Integer>::type
{
#ifdef T8_WITH_NETCDF
  T8_ASSERT (exponent >= 0);
  /* Declare the result variable */
  Integer result = 1;

  /* Perform exponentiation by squaring */
  while (true) {
    /* Check whether the (current) exponent is odd */
    if (exponent & 1) {
      /* Multiply the base */
      result *= base;
    }
    /* Halve the exponent */
    exponent >>= 1;
    /* Check whether the exponentiation continues */
    if (!exponent) {
      break;
    }
    /* Square the base */
    base *= base;
  }
  /* Return the result of the integer exponentiation */
  return result;
#endif
}

#endif /* !T8_NC_UTILITIES_HXX */
