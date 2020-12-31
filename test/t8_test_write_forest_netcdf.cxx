#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
#include <netcdf_par.h>
#endif
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_netcdf.h>
#include <t8_netcdf.h>

T8_EXTERN_C_BEGIN ();

void
t8_example_netcdf_write_forest (sc_MPI_Comm comm, int mpirank)
{
#if T8_WITH_NETCDF
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_gloidx_t         num_elements;
  int                 level = 0;

  /* Create a default scheme */
  default_scheme = t8_scheme_new_default_cxx ();

  /* Construct a cube coarse mesh */
  //cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  /*Construct a 3D hybrid hypercube as a cmesh */
  cmesh = t8_cmesh_new_hypercube_hybrid (3, comm, 0, 0);

  t8_global_productionf ("New cmesh was created\n");

  /* Build a (partioined) uniform forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);

  t8_global_productionf ("New forest was created\n");

  /* Print out the number of local elements of each  process */
  num_elements = t8_forest_get_num_element (forest);
  printf ("Rank %d has %ld elements\n", mpirank, num_elements);

  /* *Example user-defined NetCDF variable* */
  /* Allocate the data which lays on the several processes */
  /* Those user-defined variables are currently only meant to maintain a single value per (process-local) element */
  int                *var_rank = (int *) T8_ALLOC (int, num_elements);
  for (int j = 0; j < num_elements; j++) {
    var_rank[j] = mpirank;
  }
  /* Create an extern INT -NetCDF variable and receive the pointer to it */
  /* Currently INT-NetCDF and DOUBLE-NetCDF variables are possible */
  t8_netcdf_variable_t *ext_var_mpirank =
    t8_netcdf_variable_int_init ("mpirank",
                                 "Mpirank which the element lays on",
                                 "integer", var_rank);

  /* Create an array of pointers to extern NetCDF-variables, further extern NetCDF-Variables could be created and appended to the array */
  t8_netcdf_variable_t **ext_vars = new t8_netcdf_variable_t *[1];
  ext_vars[0] = ext_var_mpirank;

  /* Write the forest to NetCDF */
  t8_forest_write_netcdf (forest, "TryForestNetCDFParallelWithExtVar",
                          "Example Uniform Forest", 3, 1, ext_vars);

  t8_global_productionf ("The forest has been written to a NetCDF file\n");

  /* Destroy the forest */
  t8_forest_unref (&forest);

  /* Free the allocated array of pointers to extern NetCDF-variables */
  delete[]ext_vars;

  /* Free the allocated memory of the extern NetCDF-variables which was created by calling the 'destroy' function */
  t8_netcdf_variable_destroy (ext_var_mpirank);

  /* Free the data of the user-defined variable */
  T8_FREE (var_rank);

#endif
}

int
main (int argc, char **argv)
{
#if T8_WITH_NETCDF
  int                 mpiret, mpisize, mpirank;

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

  /* Call to an example function which writes out a forest in NetCDF-Format */
  t8_example_netcdf_write_forest (sc_MPI_COMM_WORLD, mpirank);

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

#endif
  return 0;
}

T8_EXTERN_C_END ();
