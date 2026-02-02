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

/** \file t8_netcdf.h
 * This Header-File holds two typedef integer datatypes for the data of NetCDF integer variables, as well as a struct and functions to create extern elementwise data-variables which should be written additionally to a Forest or Cmesh NetCDF-File.
 */

#ifndef T8_NETCDF_H
#define T8_NETCDF_H

#include <t8.h>

T8_EXTERN_C_BEGIN ();

/** Datatype for 64 bit NetCDF Data */
typedef int64_t t8_nc_int64_t;

/** Datatype for 32 bit NetCDF Data */
typedef int32_t t8_nc_int32_t;

/** This enumeration contains all possible netCDF variable datatypes (int, int64, double). */
typedef enum t8_netcdf_variable_type {
  /** Symbolizes netCDF variable datatype which holds 32-bit integer data */
  T8_NETCDF_INT = 0,
  /** Symbolizes netCDF variable datatype which holds 64-bit integer data */
  T8_NETCDF_INT64 = 1,
  /** Symbolizes netCDF variable datatype which holds double data */
  T8_NETCDF_DOUBLE = 2
} t8_netcdf_variable_type_t;

/** Struct for elementwise data variable for a NetCDF file */
typedef struct
{
  /** short name of the variable*/
  const char *variable_name;
  /** long name of the variable*/
  const char *variable_long_name;
  /** unit of the variable*/
  const char *variable_units;
  /** datatype of the variable*/
  t8_netcdf_variable_type_t datatype;
  /** The unique identifier netCDF assigns to this variable once it is defined in the file*/
  int var_user_dimid;
  /** user data for the variable*/
  sc_array_t *var_user_data;
} t8_netcdf_variable_t;

/** Create an extern double variable which additionally should be put out to the NetCDF File
 * \param [in]  var_type    Defines the datatype of the variable, either T8_NETCDF_INT, T8_NETCDF_INT64 or T8_NETCDF_DOUBLE.
 * \param [in]  var_name    A String which will be the name of the created variable.
 * \param [in]  var_long_name    A string describing the variable a bit more and what it is about.
 * \param [in]  var_unit    The units in which the data is provided.
 * \param [in]  var_data    A sc_array_t holding the elementwise data of the variable.
 */
t8_netcdf_variable_t *
t8_netcdf_create_var (t8_netcdf_variable_type_t var_type, const char *var_name, const char *var_long_name,
                      const char *var_unit, sc_array_t *var_data);

/** Create an extern integer variable which additionally should be put out to the NetCDF File (The distinction if it will be a NC_INT or NC_INT64 variable is based on the elementsize of the given sc_array_t)
 * \param [in]  var_name    A String which will be the name of the created variable.
 * \param [in]  var_long_name    A string describing the variable a bit more and what it is about.
 * \param [in]  var_unit    The units in which the data is provided.
 * \param [in]  var_data    A sc_array_t holding the elementwise data of the variable.
 */
t8_netcdf_variable_t *
t8_netcdf_create_integer_var (const char *var_name, const char *var_long_name, const char *var_unit,
                              sc_array_t *var_data);

/** Create an extern double variable which additionally should be put out to the NetCDF File
 * \param [in]  var_name    A String which will be the name of the created variable.
 * \param [in]  var_long_name    A string describing the variable a bit more and what it is about.
 * \param [in]  var_unit    The units in which the data is provided.
 * \param [in]  var_data    A sc_array_t holding the elementwise data of the variable.
 */
t8_netcdf_variable_t *
t8_netcdf_create_double_var (const char *var_name, const char *var_long_name, const char *var_unit,
                             sc_array_t *var_data);

/** Free the allocated memory of the a t8_netcdf_variable_t
 * \param [in]  var_destroy    A t8_netcdf_t variable whose allocated memory should be freed.
 */
/* Free the allocated NetCDF variable */
void
t8_netcdf_variable_destroy (t8_netcdf_variable_t *var_destroy);

T8_EXTERN_C_END ();

#endif /* T8_NETCDF_H */
