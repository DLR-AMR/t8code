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

/** \file t8_forest_netcdf.h
 * The Header-File providing a function to write out a forest in a NetCDF-4 file.
 */

#ifndef T8_FOREST_NETCDF_H
#define T8_FOREST_NETCDF_H

#include <t8_forest.h>
#include <t8_netcdf.h>

T8_EXTERN_C_BEGIN ();

/** Creates a netCDF-4 file containing the (geometrical) information about the given forest mesh and additional elementwise data variables
 * This functions calls 't8_forest_write_netcdf_ext()' with default values for 'netcdf_var_storage_mode' (NC_CONTIGUOUS), 'netcdf_var_mpi_access' (NC_INDEPENDENT), 'coordinate_chunksize' (nullptr), fill_mode (NC_NOFILL), cmode (NC_NETCDF4), and multifile_mode (false).
 * \param [in]  forest    A forest.
 * \param [in]  file_prefix    A string which holds the file's name (output file will be 'file_prefix.nc').
 * \param [in]  file_title    A string to caption the NetCDF-File.
 * \param [in]  dim    The Dimension of the forest mesh (2D or 3D).
 * \param [in]  num_extern_netcdf_vars    The number of extern user-defined variables which hold elementwise data (if none, set it to 0).
 * \param [in]  ext_variables An array of pointers of the herein before mentioned user-defined variables (if none, set it to NULL).
 * \param [in]  comm The sc_MPI_Communicator to use.
 */
void                t8_forest_write_netcdf (t8_forest_t forest,
                                            const char *file_prefix,
                                            const char *file_title, int dim,
                                            int num_extern_netcdf_vars,
                                            t8_netcdf_variable_t *
                                            ext_variables[],
                                            sc_MPI_Comm comm);

/** Creates a netCDF-4 file containing the (geometrical) information about the given forest mesh and additional elementwise data variables
 * \param [in]  forest    A forest.
 * \param [in]  file_prefix    A string which holds the file's name (output file will be 'file_prefix.nc').
 * \param [in]  file_title    A string to caption the NetCDF-File.
 * \param [in]  dim    The Dimension of the forest mesh (2D or 3D).
 * \param [in]  num_extern_netcdf_vars    The number of extern user-defined variables which hold elementwise data (if none, set it to 0).
 * \param [in]  ext_variables An array of pointers of the herein before mentioned user-defined variables (if none, set it to NULL).
 * \param [in]  comm The sc_MPI_Communicator to use.
 * \param [in]  netcdf_var_storage_mode This defines whether the variables in the netCDF-4 File are stored contiguously (NC_CONTIGUOUS) or in chunks (NC_CHUNKED).
 * \param [in]  coordinate_chunksize a pointer to a single size_t indicating the chunksize used only for the coordinate variables. Ignored if netcdf_var_storage_mode is not NC_CHUNKED. May be nullptr indicating the netcdf default
 * \param [in]  netcdf_var_mpi_access This defines whether the mpi processes are writing the variables' data collectively or independently.
 * \param [in]  fill_mode the fillmode used for the netcdf file. Must be either NC_NOFILL or NC_FILL
 * \param [in]  cmode the given flags are potentially combined with other flags (eg. NC_CLOBBER) and given to nc_create[_par]. Supported modes are NC_NETCDF4 aka netCDF4/HDF5 and NC_64BIT_DATA aka CDF-5.
 * \param [in]  multifile_mode if true, every MPI rank writes its local data to its own file, which can have performance benefits but leaves you with many files. In this mode, netcdf_var_mpi_access is irrelevant and ignored.
 * \note Only use this extended version of \ref t8_forest_write_netcdf if you know exactly what you are doing. In all other cases, we recommend using \ref t8_forest_write_netcdf instead.
 */
void                t8_forest_write_netcdf_ext (t8_forest_t forest,
                                                const char *file_prefix,
                                                const char *file_title,
                                                int dim,
                                                int num_extern_netcdf_vars,
                                                t8_netcdf_variable_t *
                                                ext_variables[],
                                                sc_MPI_Comm comm,
                                                int netcdf_var_storage_mode,
                                                const size_t
                                                *coordinate_chunksize,
                                                int netcdf_var_mpi_access,
                                                int fill_mode, int cmode,
                                                bool multifile_mode);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_NETCDF_H */
