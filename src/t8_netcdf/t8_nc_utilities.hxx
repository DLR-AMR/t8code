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

#include <t8.h>
#include <type_traits>

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
