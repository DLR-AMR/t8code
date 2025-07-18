#include <t8_mra/t8_mra.hpp>
#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_geometry/t8_geometry_with_vertices.h"

t8_cmesh_t
t8_cmesh_new_debugging (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[18] = {
    0, 0, 0, 0, 1, 0, 1, 0, 0,  //triangle 1
    1, 1, 0, 1, 0, 0, 0, 1, 0,  //triangle 2
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

int
main (int argc, char** argv)
{
  int mpiret;
  sc_MPI_Comm comm;
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);
  comm = sc_MPI_COMM_WORLD;

  printf ("Init done\n");

  int max_level = 8;
  double c_thresh = 1.0;
  int dunavant_rule = 10;

  constexpr int P = 3;
  constexpr int U = 1;
  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;

  printf ("\nSettings:\n max_level %d\nc_thresh %f\ndunavant rule %d\nP %d\nU %d\nelem_type %d\n", max_level, c_thresh,
          dunavant_rule, P, U, T8_ECLASS_TRIANGLE);

  auto* test_scheme = t8_scheme_new_default ();
  t8_cmesh_t cmesh = t8_cmesh_new_debugging (comm);
  printf ("created test_scheme and cmesh\n");

  t8_mra::levelindex_map<element_data_type>* lmi_map = new t8_mra::levelindex_map<element_data_type> ();
  printf ("created lmi_map\n");

  t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra_test (max_level, c_thresh, dunavant_rule, comm);
  printf ("created mra object\n");

  /// Velis debugging example
  auto f4 = [] (double x, double y) {
    //if ((x == -1.) && (y == -1.)) return 6.;
    if (x < 0.41)
      return 0.;
    double r4 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
    double r = sqrt (r4);
    if (r > 1. / 3.)
      return 0.;
    r *= 3.;
    r4 *= 9.;
    r4 *= r4;
    double rm1 = r - 1.;
    double rm1h2 = rm1 * rm1;
    double rm1h3 = rm1 * rm1h2;
    return 1. - r4 + 4. * r4 * rm1 - 10. * r4 * rm1h2 + 20 * r4 * rm1h3;
  };

  // auto test_forest
  //   = mra_test.initialize_data (lmi_map, cmesh, test_scheme, 1u, [] (double x, double y) { return x + y; });
  auto forest = mra_test.initialize_data (lmi_map, cmesh, test_scheme, 1u, f4);
  printf ("initialize data\n");

  auto* user_data = reinterpret_cast<t8_mra::forest_data<element_data_type>*> (t8_forest_get_user_data (forest));

  printf ("size init data: %zu\n", user_data->lmi_map->size ());

  sc_array_destroy (user_data->lmi_idx);
  delete lmi_map;

  t8_forest_unref (&forest);
  T8_FREE (user_data);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  printf ("freed everything...\n");

  return 0;
}
