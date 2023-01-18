#include <t8_cmesh/t8_cmesh_examples.h>

int
main (int argc, char **argv)
{
  t8_cmesh_t          cmesh;
  p4est_connectivity_t *conn;

  SC_CHECK_MPI(sc_MPI_Init (&argc, &argv));

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);
  
  conn = p4est_connectivity_new_brick (1, 1, 0, 0);
  cmesh = t8_cmesh_new_from_p4est (conn, sc_MPI_COMM_WORLD, 0);
  p4est_connectivity_destroy (conn);
  t8_cmesh_destroy (&cmesh);

  sc_finalize ();
  SC_CHECK_MPI (sc_MPI_Finalize ());
  return 0;
}
