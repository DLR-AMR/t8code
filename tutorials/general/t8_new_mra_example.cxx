#include <t8_mra/t8_mra.hpp>

#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_geometry/t8_geometry_with_vertices.h"
#include "t8_mra/data/cell_data.hpp"
#include "t8_vtk.h"

t8_cmesh_t
t8_cmesh_square (sc_MPI_Comm comm)
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

t8_cmesh_t
t8_cmesh_l_shape (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[36] = {
    0.5, 0.5, 0, 1, 0, 0, 1,   0.5, 0,  //triangle 1
    0.5, 0.5, 0, 1, 0, 0, 0,   0,   0,  //triangle 2
    0.5, 0.5, 0, 0, 1, 0, 0,   0,   0,  //triangle 3
    0.5, 0.5, 0, 0, 1, 0, 0.5, 1,   0,  //triangle 4
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 2, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

template <typename T>
void
t8_write_vtu (t8_forest_t forest, const char* prefix)
{
  const auto total_num_elements = t8_forest_get_global_num_leaf_elements (forest);

  constexpr auto U_DIM = T::U_DIM;
  auto num_fields_to_write = U_DIM;

  std::array<double*, U_DIM> element_data;
  for (auto k = 0u; k < U_DIM; ++k)
    element_data[k] = T8_ALLOC (double, total_num_elements);

  std::array<t8_vtk_data_field_t, 3> vtk_data;
  for (auto k = 0u; k < U_DIM; ++k) {
    vtk_data[k].type = T8_VTK_SCALAR;
    strcpy (vtk_data[k].description, ("u" + std::to_string (k)).c_str ());
    vtk_data[k].data = element_data[k];
  }

  const t8_element_t* element;
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  auto current_index = 0u;
  for (auto tree_idx = 0u, current_index = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    auto* data = t8_mra::get_mra_forest_data<T> (forest);
    for (auto ele_idx = 0u; ele_idx < num_elements; ++ele_idx, ++current_index) {
      element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);

      const auto lmi = t8_mra::get_lmi_from_forest_data<T> (data, current_index);
      const auto mean_val = t8_mra::mean_val<T> (forest, tree_idx, lmi, element);

      for (auto k = 0u; k < U_DIM; ++k)
        element_data[k][current_index] = mean_val[k];
    }
  }

  int write_treeid = 1;
  int write_mpirank = 0;
  int write_level = 1;
  int write_element_id = 1;
  int write_ghosts = 0;
  t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts, 0,
                           0, num_fields_to_write, vtk_data.data ());

  for (auto k = 0u; k < U_DIM; ++k)
    T8_FREE (element_data[k]);
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

  /// Velis example

  /// Half circle + jump middle
  auto f4 = [] (double x, double y) -> std::array<double, 2> {
    if (x < 0.41)
      return { 0.0, 0.0 };
    auto r4 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
    auto r = std::sqrt (r4);

    if (r > 1.0 / 3.0)
      return { 0.0, 0.0 };

    r *= 3.0;
    r4 *= 9.0;
    r4 *= r4;

    const auto rm1 = r - 1.0;
    const auto rm1h2 = rm1 * rm1;
    const auto rm1h3 = rm1 * rm1h2;

    return { 1.0 - r4 + 4.0 * r4 * rm1 - 10.0 * r4 * rm1h2 + 20.0 * r4 * rm1h3,
             -(1.0 - r4 + 4.0 * r4 * rm1 - 10.0 * r4 * rm1h2 + 20.0 * r4 * rm1h3) };
  };

  /// Different oscillations
  auto f5 = [] (double x, double y) -> std::array<double, 1> { return { std::sin (1 / (1.001 - x * y)) }; };

  /// Quarter circle + jump along circle
  auto f3 = [] (double x, double y) -> std::array<double, 1> {
    double r = x * x + y * y;
    return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
  };

  auto max_level = 7u;
  auto c_thresh = 0.8;
  auto gamma = 1.0;  /// Order of convergence
  auto dunavant_rule = 10;

  bool balanced = true;

  constexpr int P = 4;
  constexpr int U = 1;

  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;
  using mra_type = t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P>;

  auto* test_scheme = t8_scheme_new_default ();
  // t8_cmesh_t cmesh = t8_cmesh_square (comm);
  t8_cmesh_t cmesh = t8_cmesh_l_shape (comm);
  printf ("created test_scheme and cmesh\n");

  mra_type mra_test (max_level, c_thresh, gamma, dunavant_rule, balanced, comm);
  printf ("created mra object\n");

  mra_test.initialize_data (cmesh, test_scheme, max_level, f3);
  printf ("Initialize data\n");
  printf ("Size init data: %zu\n", mra_test.get_lmi_map ()->size ());

  // t8_write_vtu<element_data_type> (mra_test.forest, ("reference_f3_" + std::to_string (init_level)).c_str ());

  mra_test.coarsening (0, max_level);

  printf ("Size after coarsening: %zu\n", mra_test.get_lmi_map ()->size ());

  t8_write_vtu<element_data_type> (
    mra_test.forest,
    ("f3_l_shape_balanced_coarsening_P" + std::to_string (P) + "_" + std::to_string (max_level)).c_str ());

  printf ("Finished writing file\n");

  mra_test.cleanup ();
  printf ("cleaned everything...\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  printf ("Program finished...\n");
  return 0;
}
