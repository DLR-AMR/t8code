#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_vtk/t8_vtk_writer.h>

int
main (int argc, char **argv)
{
  int mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  const char *fileprefix = "examplefiles/mwe";
  const char *outfileprefix = "t8_mwe";

  t8_cmesh_t cmesh = t8_cmesh_from_msh_file(fileprefix, 0, sc_MPI_COMM_WORLD, 2, 0, 0);
  t8_cmesh_vtk_write_file(cmesh, outfileprefix);
  t8_cmesh_destroy(&cmesh);

  const char *fileprefix_rev = "examplefiles/mwe_reversed";
  const char *outfileprefix_rev = "t8_mwe_reversed";

  t8_cmesh_t cmesh_rev = t8_cmesh_from_msh_file(fileprefix_rev, 0, sc_MPI_COMM_WORLD, 2, 0, 0);
  t8_cmesh_vtk_write_file(cmesh_rev, outfileprefix_rev);
  t8_cmesh_destroy(&cmesh_rev);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
