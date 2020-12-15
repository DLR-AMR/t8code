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
#include <netcdf_par.h>
/* Standard netcdf error function */
#define ERRCODE 2
#define ERR(e) {t8_global_productionf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#endif
#include <t8_cmesh.h>
#include <t8_cmesh_netcdf.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_eclass.h>

typedef struct
{
  int                 ncidp;
  int                 dim_id;
  int                 elem_id;
} just_data_t;

void
t8_cmesh_write_netcdf_parallel_try (t8_cmesh_t cmesh, sc_MPI_Comm comm,
                                    just_data_t * context, int rank)
{
  //int ncidp;
  int                 retval;
  //int nMesh_elem_dimid;
  int                 nMesh_elem = 6;
  const char         *filename = "try_parallelJetzt.nc";

  if ((retval =
       nc_create_par (filename, NC_CLOBBER | NC_NETCDF4, comm,
                      sc_MPI_INFO_NULL, &context->ncidp))) {
    ERR (retval);
  }
  if ((retval =
       nc_def_dim (context->ncidp, "number_elements", nMesh_elem,
                   &context->dim_id))) {
    ERR (retval);
  }
  if ((retval =
       nc_def_var (context->ncidp, "elements",
                   NC_DOUBLE, 1, &context->dim_id, &context->elem_id))) {
    ERR (retval);
  }
#if 0
  if ((retval =
       nc_var_par_access (context->ncidp, context->elem_id, NC_COLLECTIVE))) {
    ERR (retval);
  }
#endif
  if ((retval = nc_enddef (context->ncidp))) {
    ERR (retval);
  }
  if (rank == 0) {
    int                 daten[3] = { 1, 2, 3 };
    size_t              startp = 0;
    size_t              countp = 3;
    if ((retval =
         nc_put_vara_int (context->ncidp, context->elem_id, &startp, &countp,
                          &daten[0]))) {
      ERR (retval);
    }
  }
  else if (rank == 2) {
    int                 daten2[3] = { 4, 5, 6 };
    size_t              startp2 = 3;
    size_t              countp2 = 2;
    if ((retval =
         nc_put_vara_int (context->ncidp, context->elem_id, &startp2,
                          &countp2, &daten2[0]))) {
      ERR (retval);
    }
  }
  if ((retval = nc_close (context->ncidp))) {
    ERR (retval);
  }

}

int
main (int argc, char **argv)
{
  int                 mpiret, mpisize, mpirank;
  t8_cmesh_t          cmesh;

#if 1
  /* Initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Initialize sc */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code */
  t8_init (SC_LP_PRODUCTION);

  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
  SC_CHECK_MPI (mpiret);
#endif

#if 0
  if (mpirank == 0) {
#if 0
    char               *mesh_name2 = "NewTry2DVersion2";
    cmesh =
      t8_cmesh_from_msh_file ("test/test_msh_file_vers2_ascii", 0,
                              sc_MPI_COMM_WORLD, 2, 0);
    t8_cmesh_write_netcdf (cmesh, mesh_name2, "BeispielDaten", 2);
    retval = t8_cmesh_write_netcdf2D (cmesh, mesh_name2, "BeispielDaten");
#endif
#if 0
    //try_netcdf_other_format();
    cmesh =
      t8_cmesh_from_msh_file ("test/test_msh_file_vers2_ascii", 0,
                              sc_MPI_COMM_WORLD, 2, 0);
    char               *mesh_name = "NewTry2DWithChangesJetzt";
    t8_cmesh_write_netcdf (cmesh, mesh_name, "Beispiel 2D", 2);
    t8_global_productionf ("NetCDF output\n");
    t8_cmesh_destroy (&cmesh);
#endif

  }
#endif

#if 1
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_WORLD, 0, 1, 0);

  printf ("Number of local trees on process %d : %d\n", mpirank,
          (int) t8_cmesh_get_num_local_trees (cmesh));
//t8_cmesh_write_netcdf_parallel_try(cmesh, sc_MPI_COMM_WORLD, &content, mpirank);
  char               *mesh_name = "NewTry3DParallel";
  t8_cmesh_write_netcdf (cmesh, mesh_name, "Beispiel 3D Parallel", 3);
  t8_global_productionf ("NetCDF output\n");
  t8_cmesh_destroy (&cmesh);
#endif
#if 1
  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
#endif
  return 0;
}
