#include <gtest/gtest.h>
#include <t8.h>

int
main (int argc, char **argv)
{

  int mpiret;
  sc_MPI_Comm mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (SC_LP_DEFAULT);

  ::testing::InitGoogleTest (&argc, argv);

  int retval = RUN_ALL_TESTS ();

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return retval;
}
