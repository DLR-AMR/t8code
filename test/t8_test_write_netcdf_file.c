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
/* Copied from the linkage netcdf test. */
/* In this test we create a netcdf in memory file and close it.
 * The purpose of this test is to check whether t8code successfully links
 * against netcdf.
 * If t8code was not configured with --with-netcdf then this test
 * does nothing and is always passed.
 */

#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) {t8_global_productionf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif
#include <t8_cmesh.h>
#include <t8_cmesh_netcdf.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_eclass.h>
#include <t8_netcdf.h>

void
t8_example_netcdf_write_cmesh (int mpirank)
{

  t8_cmesh_t          cmesh;

  /* Construct a (partitioned) cmesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_WORLD, 0, 1, 0);
  //cmesh = t8_cmesh_new_hypercube_hybrid(3, sc_MPI_COMM_WORLD, 1, 0);

  /* Number of process local elements */
  printf ("Number of local trees on process %d : %d\n", mpirank,
          t8_cmesh_get_num_local_trees (cmesh));

  /* *Example user-defined NetCDF variable* */
  /* Allocate the data which lays on the several processes */
  /* Those user-defined variables are currently only meant to maintain a single value per (process-local) element */
  int                *var_rank =
    (int *) T8_ALLOC (int, t8_cmesh_get_num_local_trees (cmesh));
  for (int j = 0; j < t8_cmesh_get_num_local_trees (cmesh); j++) {
    var_rank[j] = mpirank;
  }
  /* Create an extern INT -NetCDF variable and receive the pointer to it */
  /* Currently INT-NetCDF and DOUBLE-NetCDF variables are possible */
  t8_netcdf_variable_t *ext_var_mpirank =
    t8_netcdf_variable_int_init ("mpirank",
                                 "Mpirank which the element lays on",
                                 "integer", var_rank);

  /* Create an array of pointers to extern NetCDF-variables, further extern NetCDF-Variables could be created and appended to the array */
  t8_netcdf_variable_t *ext_vars[1] = { ext_var_mpirank };

  char               *mesh_name = "NewCmesh3DParallel";

  /* Write the cmesh to NetCDF */
  t8_cmesh_write_netcdf (cmesh, mesh_name, "Beispiel 3D Parallel Cmesh", 3, 1,
                         ext_vars);
  t8_global_productionf ("NetCDF output of the cmesh has been written\n");

  /* Destroy the cmesh */
  t8_cmesh_destroy (&cmesh);

  /* Free the allocated memory of the extern NetCDF-variables which was created by calling the 'destroy' function */
  t8_netcdf_variable_destroy (ext_var_mpirank);

  /* Free the data of the user-defined variable */
  T8_FREE (var_rank);

}

int
main (int argc, char **argv)
{
#if T8_WITH_NETCDF
  int                 mpiret, mpirank;

  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);

  /* An example function that outputs a cmesh in NetCDF-Format */
  t8_example_netcdf_write_cmesh (mpirank);

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
#endif
  return 0;
}
