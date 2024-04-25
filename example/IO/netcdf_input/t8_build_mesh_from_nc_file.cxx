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

#include <t8_netcdf/t8_nc.hxx>
#include <t8_netcdf/t8_nc_utilities.hxx>
#include <t8_netcdf/t8_nc_hyperslab.hxx>
#include <t8_netcdf/t8_nc_mesh.hxx>

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  t8_global_essentialf ("NetCDF Input Example\n");

  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 96);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 64);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 1);

  /* Make a hyperslab from all dimensions */
  const t8_nc_hyperslab_t hyperslab (lon_dimension, lat_dimension, lev_dimension);

  /* Open a netCDF file for reading */
  t8_nc_data_t nc_data ("../../data/test_nc4_file.nc", t8_nc_opening_mode_t::SERIAL);

  /* Inquire the data of the considered variables */
  nc_data.inquire_variables (hyperslab, "p2t");

  /* Close the file after the data has been read */
  nc_data.close_file_handle ();

  /* Seize the inquired data */
  std::vector<InputVar> nc_variables = nc_data.transfer_data ();
  t8_productionf ("Size of InputVar vector: %ld\n", nc_variables.size ());

  //auto [initial_embedded_forest, current_max_refinement_lvl_initial] = t8_nc_build_initial_embedded_mesh();
  //auto [initial_congruent_forest, current_max_refinement_lvl_congruent] = t8_nc_build_initial_congruent_mesh();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
