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

#include <t8_netcdf/t8_nc.h>
#include <t8_netcdf/t8_nc_data.hxx>
int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  t8_global_essentialf ("NetCDF Input Example\n");

  t8_universal_type_t uni;

  uni = static_cast<double> (2.0);

  uni = static_cast<int32_t> (344);

  t8_nc_data_t nc_data
    = t8_nc_start ("../../data/test_nc4_file.nc", t8_nc_opening_mode::T8_NC_PARALLEL, sc_MPI_COMM_WORLD);

  const char *var_names[] = { "p2t" };
  const int num_variables = 1;

  const size_t start_ptr[3] = { 0, 0, 0 };     //Example netCDF File
  const size_t count_ptr[3] = { 1, 73, 144 };  //Example netCDF File

  const enum t8_nc_geo_mesh_type mesh_type = T8_NC_CONGRUENT_MESH;
  const enum t8_nc_geo_mesh_form mesh_form = T8_NC_RECTANGULAR;
  const enum t8_nc_geo_mesh_elements mesh_elems = T8_NC_QUAD_ELEMENTS;

  //int num_procs_per_dim[3] = {1,1,2};
  //t8_nc_set_hint_read_data_blockwise_in_parallel(nc_data, 3, num_procs_per_dim);

  t8_nc_construct_mesh_for_variables (nc_data, num_variables, var_names, start_ptr, count_ptr, mesh_type, mesh_form,
                                      mesh_elems);

  t8_nc_finish (nc_data);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
