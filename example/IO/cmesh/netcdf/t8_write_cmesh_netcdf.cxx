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
/*
  In this example, the use of functions in order to write out a cmesh in the netCDF format is exemplarily displayed.
  Additionally to the UGRID representation of the cmesh, element-wise user-defined variables are created and also written to the netCDF file.
 */

#include <t8.h>
#include <netcdf.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) \
  { \
    t8_global_productionf ("Error: %s\n", nc_strerror (e)); \
    exit (ERRCODE); \
  }
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_netcdf.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include <t8_element/t8_eclass.h>
#include <t8_netcdf/t8_netcdf.h>

void
t8_example_netcdf_write_cmesh (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;
  t8_nc_int32_t *var_rank;
  sc_array_t *var_ranks;
  sc_array_t *var_random_values;
  double *random_values;
  t8_netcdf_variable_t *ext_var_mpirank;
  t8_netcdf_variable_t *ext_var_random_values;
  int mpiret;
  int mpirank;
  int j;

  /* Receive the process-local MPI rank */
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* Construct a hybrid cmesh */
  cmesh = t8_cmesh_new_hypercube_hybrid (comm, 1, 0);

  /* Number of process-local elements */
  t8_global_productionf ("Number of local trees on process %d : %d\n", mpirank, t8_cmesh_get_num_local_trees (cmesh));

  /* *Example user-defined NetCDF variable, mpirank* */
  /* Allocate the data which lays on the several processes */

  /* Those user-defined variables are currently only meant to maintain a single value per (process-local) element */
  var_rank = T8_ALLOC (t8_nc_int32_t, t8_cmesh_get_num_local_trees (cmesh));
  /* Write out the mpirank of each (process-local) element */
  for (j = 0; j < t8_cmesh_get_num_local_trees (cmesh); j++) {
    var_rank[j] = mpirank;
  }

  /* Create a new sc_array_t which provides the data for the NetCDF variables, in this case the mpirank each element lays on */
  var_ranks = sc_array_new_data (var_rank, sizeof (t8_nc_int32_t), t8_cmesh_get_num_local_trees (cmesh));

  /* Create an integer NetCDF variable; parameters are (name of the variable, descriptive long name of the variable, description of the data's unit, pointer to sc_array_t which provides the data) */
  ext_var_mpirank = t8_netcdf_create_integer_var ("mpirank", "Mpirank which the element lays on", "integer", var_ranks);

  /* *Example user-defined NetCDF variable, random values* */
  /* Create random values */
  random_values = T8_ALLOC (double, t8_cmesh_get_num_local_trees (cmesh));

  for (j = 0; j < t8_cmesh_get_num_local_trees (cmesh); j++) {
    random_values[j] = rand () / (double) rand ();
  }
  /* Create a new sc_array_t which provides the data for the NetCDF variables, in this case just random values */
  var_random_values = sc_array_new_data (random_values, sizeof (double), t8_cmesh_get_num_local_trees (cmesh));

  /* Create a double NetCDF variable; parameters are (name of the variable, descriptive long name of the variable, description of the data's unit (i.e. degrees Celsius), pointer to sc_array_t which provides the data) */
  ext_var_random_values
    = t8_netcdf_create_double_var ("random_values", "Random values in [0,10)", "double", var_random_values);

  /* Create an array of pointers to extern NetCDF-variables, further extern NetCDF-Variables could be created and appended to the array */
  t8_netcdf_variable_t *ext_vars[2] = { ext_var_mpirank, ext_var_random_values };

  /* Name of the NetCDF-File */
  const char *mesh_name = "T8_Example_NetCDF_Cmesh";

  /* Write the cmesh to NetCDF */
  t8_cmesh_write_netcdf (cmesh, mesh_name, "Example 3D parallel cmesh", 3, 2, ext_vars, comm);

  t8_global_productionf ("NetCDF output of the cmesh has been written.\n");

  /* Destroy the cmesh */
  t8_cmesh_destroy (&cmesh);

  /* Free the allocated memory of the extern NetCDF-variables which was created by calling the 'destroy' function */
  t8_netcdf_variable_destroy (ext_var_mpirank);
  t8_netcdf_variable_destroy (ext_var_random_values);

  /* Destroy the allocated sc_array_t */
  sc_array_destroy (var_ranks);
  sc_array_destroy (var_random_values);

  /* Free the data of the user-defined variable */
  T8_FREE (var_rank);
  T8_FREE (random_values);
}

int
main (int argc, char **argv)
{
  int mpiret;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  /* An example function that outputs a cmesh in NetCDF-Format */
  t8_example_netcdf_write_cmesh (sc_MPI_COMM_WORLD);

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
