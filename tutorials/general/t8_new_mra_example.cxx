#include <t8_mra/t8_mra.hpp>
#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_forest/t8_forest_geometrical.h"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_geometry/t8_geometry_with_vertices.h"
#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_vtk.h"

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

template <typename T>
void
t8_write_vtu (t8_forest_t forest, t8_mra::forest_data<T>* data, const char* prefix)
{
  const auto total_num_elements = t8_forest_get_global_num_leaf_elements (forest);
  double* element_data = T8_ALLOC (double, total_num_elements);

  auto num_data = 1;
  t8_vtk_data_field_t vtk_data;
  vtk_data.type = T8_VTK_SCALAR;
  strcpy (vtk_data.description, "Element own data");
  vtk_data.data = element_data;

  /// TODO Easy access to element in forest
  auto get_value = [&] (const t8_mra::forest_data<T>* forest_data, auto idx) {
    return *((t8_mra::levelmultiindex<T::Shape>*) t8_sc_array_index_locidx ((forest_data->lmi_idx), idx));
  };

  const t8_element_t* element;
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  auto current_index = 0u;
  for (auto tree_idx = 0u, current_index = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
      const auto vol = t8_forest_element_volume (forest, tree_idx, element);

      const auto lmi = get_value (data, current_index);
      element_data[current_index] = data->lmi_map->get (lmi).u_coeffs[0];

      /// TODO Eval function
      element_data[current_index] *= t8_mra::skalierungsfunktion (0, 0.0, 0.0) * std::sqrt (1.0 / (2.0 * vol));
    }
  }

  int write_treeid = 1;
  int write_mpirank = 1;
  int write_level = 1;
  int write_element_id = 1;
  int write_ghosts = 0;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, 0,
                           0, num_data, &vtk_data);
  T8_FREE (element_data);
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

  auto f = [] (double x, double y) { return x + y; };

  printf ("Init done\n");

  auto max_level = 8u;
  auto c_thresh = 1.0;
  auto dunavant_rule = 10;

  constexpr int P = 3;
  constexpr int U = 1;
  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;
  using mra_type = t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P>;

  printf ("\nSettings:\nmax_level %d\nc_thresh %f\ndunavant rule %d\nP %d\nU %d\nelem_type %d\n", max_level, c_thresh,
          dunavant_rule, P, U, T8_ECLASS_TRIANGLE);

  auto* test_scheme = t8_scheme_new_default ();
  t8_cmesh_t cmesh = t8_cmesh_new_debugging (comm);
  printf ("created test_scheme and cmesh\n");

  mra_type mra_test (max_level, c_thresh, dunavant_rule, comm);
  printf ("created mra object\n");

  mra_test.initialize_data (cmesh, test_scheme, 8u, f4);
  // mra_test.initialize_data (cmesh, test_scheme, 4u, f);
  printf ("initialize data\n");

  printf ("size init data: %zu\n", mra_test.get_lmi_map ()->size ());
  t8_write_vtu<element_data_type> (mra_test.forest, mra_test.get_user_data (), "testi_test_8");

  printf ("Finished writing file\n");

  mra_test.cleanup ();
  printf ("freed everything...\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  printf ("Program finished...\n");
  return 0;
}
