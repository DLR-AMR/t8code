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

/** file t8_netcdf.h
 * File description
 */

#include <t8_netcdf.h>

t8_netcdf_variable_t *
t8_netcdf_variable_int_init (const char *var_name, const char *var_long_name,
                             const char *var_unit, int var_data[])
{
  t8_netcdf_variable_t *netcdf_variable = T8_ALLOC (t8_netcdf_variable_t, 1);
  netcdf_variable->variable_name = var_name;
  netcdf_variable->variable_long_name = var_long_name;
  netcdf_variable->datatype = 0;
  netcdf_variable->variable_units = var_unit;
  netcdf_variable->netcdf_data_int = var_data;
  return netcdf_variable;
}

t8_netcdf_variable_t *
t8_netcdf_variable_double_init (const char *var_name,
                                const char *var_long_name,
                                const char *var_unit, double var_data[])
{
  t8_netcdf_variable_t *netcdf_variable = T8_ALLOC (t8_netcdf_variable_t, 1);
  netcdf_variable->variable_name = var_name;
  netcdf_variable->variable_long_name = var_long_name;
  netcdf_variable->datatype = 1;
  netcdf_variable->variable_units = var_unit;
  netcdf_variable->netcdf_data_double = var_data;
  return netcdf_variable;
}

void
t8_netcdf_variable_destroy (t8_netcdf_variable_t * var_destroy)
{
  T8_ASSERT (var_destroy != NULL);
  T8_FREE (var_destroy);
}
