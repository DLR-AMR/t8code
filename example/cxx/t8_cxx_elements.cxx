#include <t8.h>
#include <t8_element.h>
#include <t8_default.h>
#include <t8_element_cxx.hxx>
#include <t8_default_cxx.hxx>


void
t8_cxx_elements_test_scheme_cxx ()
{
  t8_scheme_cxx_t  *cxx_default_scheme = t8_scheme_new_default_cxx ();
  t8_eclass_scheme_c *quad_scheme;

//  quad_scheme = cxx_default_scheme->eclass_schemes[T8_ECLASS_QUAD];

  t8_scheme_cxx_unref (&cxx_default_scheme);
}

int main (int argc, char * argv[])
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  t8_cxx_elements_test_scheme_cxx ();

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
