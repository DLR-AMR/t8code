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
  t8_nc_int32_t      *var_rank;
  double             *random_values;
  sc_array_t         *var_ranks;
  sc_array_t         *var_random_values;
  t8_netcdf_variable_t *ext_var_mpirank;
  t8_netcdf_variable_t *ext_var_random_values;
  int                 level = 0;
  int                 j;

  /* Create a default scheme */
  default_scheme = t8_scheme_new_default_cxx ();

  /* Construct a cube coarse mesh */
  /* Construct a 3D hybrid hypercube as a cmesh */
  cmesh = t8_cmesh_new_hypercube_hybrid (3, comm, 0, 0);

  t8_global_productionf ("New cmesh was created\n");

  /* Build a (partioined) uniform forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);

  t8_global_productionf ("New forest was created\n");

  /* Print out the number of local elements of each  process */
  num_elements = t8_forest_get_num_element (forest);
  t8debugf ("[t8] Rank %d has %ld elements\n", mpirank, num_elements);

  /* *Example user-defined NetCDF variable* */
  /* Currently, integer and double NetCDF variables are possible */

  /* Allocate the data which lays on the several processes */
  /* Those user-defined variables are currently only meant to maintain a single value per (process-local) element */
  var_rank = T8_ALLOC (t8_nc_int32_t, num_elements);
  /* Write out the mpirank of each (process-local) element */
  for (j = 0; j < num_elements; j++) {
    var_rank[j] = mpirank;
  }
  /* Create a new sc_array_t which provides the data for the NetCDF variables, in this case the mpirank each element lays on */
  var_ranks =
    sc_array_new_data (var_rank, sizeof (t8_nc_int32_t), num_elements);
  /* Create the integer NetCDF variable; parameters are (name of the variable, descriptive long name of the variable, description of the data's unit, pointer to sc_array_t which provides the data) */
  ext_var_mpirank =
    t8_netcdf_create_integer_var ("mpirank",
                                  "Mpirank which the element lays on",
                                  "integer", var_ranks);

  /* *Example user-defined NetCDF variable, random values* */
  /* Create random values */
  random_values = T8_ALLOC (double, num_elements);

  for (j = 0; j < num_elements; j++) {
    random_values[j] = rand () / (double) rand ();
  }
  /* Create a new sc_array_t which provides the data for the NetCDF variables, in this case just random values */
  var_random_values =
    sc_array_new_data (random_values, sizeof (double), num_elements);

  /* Create the integer NetCDF variable; parameters are (name of the variable, descriptive long name of the variable, description of the data's unit (i.e. degrees Celsius), pointer to sc_array_t which provides the data) */
  ext_var_random_values =
    t8_netcdf_create_double_var ("random values", "Random values in [0,10)",
                                 "double", var_random_values);

  /* Create an array of pointers to extern NetCDF-variables, further extern NetCDF-Variables could be created and appended to the array */
  t8_netcdf_variable_t **ext_vars = new t8_netcdf_variable_t *[2];
  ext_vars[0] = ext_var_mpirank;
  ext_vars[1] = ext_var_random_values;

  /* Write the forest to NetCDF */
  t8_forest_write_netcdf (forest, "TryForestNetCDFParallelWithExtVar",
                          "Example Uniform Forest", 3, 2, ext_vars,
                          sc_MPI_COMM_WORLD);

  t8_global_productionf ("The forest has been written to a NetCDF file\n");

  /* Destroy the forest */
  t8_forest_unref (&forest);

  /* Free the allocated array of pointers to extern NetCDF-variables */
  delete[]ext_vars;

  /* Free the allocated memory of the extern NetCDF-variables which was created by calling the 'destroy' function */
  t8_netcdf_variable_destroy (ext_var_mpirank);
  t8_netcdf_variable_destroy (ext_var_random_values);

  /* Destroy the allocated sc_array_t */
  sc_array_destroy (var_ranks);
  sc_array_destroy (var_random_values);

  /* Free the data of the user-defined variable */
  T8_FREE (var_rank);
  T8_FREE (random_values);

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
