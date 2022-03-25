#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>
#include <t8_forest_vtk.h>
#include <ctime>
T8_EXTERN_C_BEGIN ();



/* Removes an element of the family 
 * Requires a full uniform grid
 *   Example for level 4 Quad
 *   
 *   X X   X |   X     | X       |          
 *       X   | X   X X |   X     | X     X  
 *       X   | X   X X |   X     | X     X  
 *   X X   X |   X     | X       |          
 *   -------------------------------------
 *       X   | X   X X |   X     | X     X  
 *   X X   X |   X     | X       |          
 *   X X   X |   X     | X       |         
 *       X   | X   X X |   X     | X     X  
 *   -------------------------------------
 *   X X X X | X X X X | X X X X | X X   X  
 *   X X X X | X X X X |   X X   | X   X X  
 *   X X X X | X X X X |   X X   | X   X X  
 *   X X X X | X X X X | X X X X | X X   X  
 *   -------------------------------------
 *   X X X X | X X X X |   X X   | X   X X  
 *   X X X X | X X X X | X X X X | X X   X  
 *   X X X X | X X X X | X X X X | X X   X 
 *   X X X X | X X X X |   X X   | X   X X  
 * 
 * */
int
t8_test_adapt (t8_forest_t forest, t8_forest_t forest_from,
               t8_locidx_t which_tree, t8_locidx_t lelement_id,
               t8_eclass_scheme_c * ts, int is_family,
               int num_elements, t8_element_t * elements[])
{
  int                 level;

  //t8_debugf("[IL] %i\n", lelement_id);
  level = ts->t8_element_level (elements[0]);
  
  
  //int                 child_id, mod_id;
  //child_id = ts->t8_element_child_id (elements[0]); // 0 1 2 3
  //mod_id = lelement_id % ts->t8_element_num_children (elements[0]); // 0 1 2 3
  
  
  //ts->t8_element_num_children (elements[0]);
  srand(time(0));
  int r = rand() % 2;
  if (level > 0 && r) {
    return -2;
  }
  return 0;
}



int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;

  const int           level = 4;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  comm = sc_MPI_COMM_WORLD;

  // Build cmesh and uniform forest.
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, comm);
  T8_ASSERT (t8_forest_is_committed (forest));

  
  forest = t8_forest_new_adapt (forest, t8_test_adapt, 0, 0, NULL);
  t8_forest_write_vtk (forest, "t8_thin_out");


  t8_forest_unref (&forest);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();



