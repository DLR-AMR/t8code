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

#ifndef T8_NETCDF_H
#define T8_NETCDF_H

#include <t8.h>
//#include <t8_cmesh_netcdf.h>
//#include <t8_forest_netcdf.h>

typedef struct
{
  const char         *variable_name;
  const char         *variable_long_name;
  const char         *variable_units;
  int                 datatype;
  int                 var_user_dimid;
  union
  {
    int                *netcdf_data_int;
    double             *netcdf_data_double;
  };
} t8_netcdf_variable_t;

T8_EXTERN_C_BEGIN ();

t8_netcdf_variable_t *t8_netcdf_variable_int_init (const char *var_name,
                                                   const char *var_long_name,
                                                   const char *var_unit,
                                                   int *var_data);

t8_netcdf_variable_t *t8_netcdf_variable_double_init (const char *var_name,
                                                      const char
                                                      *var_long_name,
                                                      const char *var_unit,
                                                      double *var_data);

void                t8_netcdf_variable_destroy (t8_netcdf_variable_t *
                                                var_destroy);

T8_EXTERN_C_END ();

#endif /* T8_NETCDF_H */
