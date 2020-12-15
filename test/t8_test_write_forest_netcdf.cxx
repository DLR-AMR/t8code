#include <t8.h>
#if T8_WITH_NETCDF
#include <netcdf.h>
#endif
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_netcdf.h>

void
test_write_forest (sc_MPI_Comm comm, int level)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_scheme_cxx_t    *default_scheme;
  t8_gloidx_t         num_elements;

  default_scheme = t8_scheme_new_default_cxx ();
  /* Construct a cube coarse mesh */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_HEX, comm, 0, 0, 0);
  t8_global_productionf ("New cmesh was created\n");
  /* Build a uniform forest */
  forest = t8_forest_new_uniform (cmesh, default_scheme, level, 0, comm);
  t8_global_productionf ("New forest was created\n");
  num_elements = t8_forest_get_num_element (forest);
  t8_global_productionf ("The forest has %d elements\n", (int) num_elements);
  /* Write the forest to NetCDF */
  const char         *prefix = "TryForestNetCDF";
  const char         *title = "Example Uniform Forest";
  t8_forest_write_netcdf (forest, prefix, title, 3);
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
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

  if (mpirank == 0) {
    test_write_forest (sc_MPI_COMM_WORLD, 1);
  }

  /* Finalize sc */
  sc_finalize ();

  /* Finalize MPI */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
