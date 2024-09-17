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

#include <gtest/gtest.h>
#include <test/t8_gtest_custom_assertion.hxx>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_netcdf/t8_nc_input_variable.hxx>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <test/t8_gtest_macros.hxx>
#include <t8_netcdf/t8_nc_mesh.hxx>

//relativ einfach halten



TEST (t8_gtest_build_forest_from_geo_domain, embedded_forest)
{

  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 18);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 22);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 1);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;

  auto [initial_embedded_forest, current_max_refinement_lvl_initial] = 
              t8_nc_build_initial_embedded_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);
  
  EXPECT_EQ (t8_forest_get_global_num_elements(initial_embedded_forest), 430);

  //forest ausgeben lassen und nachschauen, ob das sinnvoll ist (wenn ungerade automatisch +1)

  t8_forest_unref(&initial_embedded_forest);


}

TEST (t8_gtest_build_forest_from_geo_domain, embedded_forest_3D)
{

  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 4);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 6);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 2);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;

  auto [initial_embedded_forest, current_max_refinement_lvl_initial] = 
              t8_nc_build_initial_embedded_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);
  
  EXPECT_EQ (t8_forest_get_global_num_elements(initial_embedded_forest), 64);

  //forest ausgeben lassen und nachschauen, ob das sinnvoll ist (wenn ungerade automatisch +1)

  t8_forest_unref(&initial_embedded_forest);


}

TEST (t8_gtest_build_forest_from_geo_domain, embedded_forest_3D_odd)
{

  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 3);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 5);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 2);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;

  auto [initial_embedded_forest, current_max_refinement_lvl_initial] = 
              t8_nc_build_initial_embedded_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);
  
  EXPECT_EQ (t8_forest_get_global_num_elements(initial_embedded_forest), 64);

  //forest ausgeben lassen und nachschauen, ob das sinnvoll ist (wenn ungerade automatisch +1)

  t8_forest_unref(&initial_embedded_forest);


}


TEST (t8_gtest_build_forest_from_geo_domain, embedded_forest_3D_quadratic)
{

  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 64);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 64);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 64);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;

  auto [initial_embedded_forest, current_max_refinement_lvl_initial] = 
              t8_nc_build_initial_embedded_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);
  
  EXPECT_EQ (t8_forest_get_global_num_elements(initial_embedded_forest), 262144);

  //forest ausgeben lassen und nachschauen, ob das sinnvoll ist (wenn ungerade automatisch +1)

  t8_forest_unref(&initial_embedded_forest);


}

TEST (t8_gtest_build_forest_from_geo_domain, congruent_mesh)
{
  
  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 18);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 22);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 1);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;
  
  auto [initial_congruent_forest, current_max_refinement_lvl_congruent] = 
              t8_nc_build_initial_congruent_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);

  EXPECT_EQ (t8_forest_get_global_num_elements(initial_congruent_forest), 396);

  EXPECT_EQ (t8_forest_get_num_global_trees(initial_congruent_forest), 99);

  t8_forest_unref(&initial_congruent_forest);

}

TEST (t8_gtest_build_forest_from_geo_domain, congruent_mesh_3D)
{
   
  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 8);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 7);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 3);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;
  
  auto [initial_congruent_forest, current_max_refinement_lvl_congruent] = 
              t8_nc_build_initial_congruent_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);

  EXPECT_EQ (t8_forest_get_global_num_elements(initial_congruent_forest), 168);

  EXPECT_EQ (t8_forest_get_num_global_trees(initial_congruent_forest), 168);

  t8_forest_unref(&initial_congruent_forest);

}


TEST (t8_gtest_build_forest_from_geo_domain, congruent_mesh_3D_power_of_two)
{
   
  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 2);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 8);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 4);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;
  
  auto [initial_congruent_forest, current_max_refinement_lvl_congruent] = 
              t8_nc_build_initial_congruent_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);

  EXPECT_EQ (t8_forest_get_global_num_elements(initial_congruent_forest), 64);

  EXPECT_EQ (t8_forest_get_num_global_trees(initial_congruent_forest), 8);

  t8_forest_unref(&initial_congruent_forest);

}

TEST (t8_gtest_build_forest_from_geo_domain, congruent_mesh_3D_quadartic)
{
   
  /* Define the ought to be considered parts of all dimensions */
  const t8_nc_dimension_interval_t lon_dimension (t8_nc_dimension_t::LON, 0, 64);
  const t8_nc_dimension_interval_t lat_dimension (t8_nc_dimension_t::LAT, 0, 64);
  const t8_nc_dimension_interval_t lev_dimension (t8_nc_dimension_t::LEV, 0, 64);

  const t8_nc_geo_domain_t global_domain(lon_dimension, lat_dimension, lev_dimension);

  const t8_nc_data_layout_t initial_layout = t8_nc_data_layout_t::LAT_LON;
  
  auto [initial_congruent_forest, current_max_refinement_lvl_congruent] = 
              t8_nc_build_initial_congruent_mesh(global_domain, initial_layout, sc_MPI_COMM_WORLD);

  EXPECT_EQ (t8_forest_get_global_num_elements(initial_congruent_forest), 262144);

  EXPECT_EQ (t8_forest_get_num_global_trees(initial_congruent_forest), 1);

  t8_forest_unref(&initial_congruent_forest);

}