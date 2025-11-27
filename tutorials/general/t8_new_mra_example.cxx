#include <t8_mra/t8_mra.hpp>
#include <t8_mra/t8_mra_vtk.hpp>

#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_geometry/t8_geometry_with_vertices.h"
#include "t8_mra/data/cell_data.hpp"
#include "t8_vtk.h"

#include <vector>
#include <cmath>

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

std::array<double, 1>
gaussian_function (double x, double y)
{
  const double x0 = 0.5;     // Center x
  const double y0 = 0.5;     // Center y
  const double sigma = 0.2;  // Standard deviation

  const double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
  const double val = std::exp (-r2 / (2.0 * sigma * sigma));

  return { val };
}

template <typename T>
void
t8_write_vtu (t8_forest_t forest, const char *prefix)
{
  const auto total_num_elements = t8_forest_get_global_num_leaf_elements (forest);

  constexpr auto U_DIM = T::U_DIM;
  auto num_fields_to_write = U_DIM;

  std::array<double *, U_DIM> element_data;
  for (auto k = 0u; k < U_DIM; ++k)
    element_data[k] = T8_ALLOC (double, total_num_elements);

  std::array<t8_vtk_data_field_t, 3> vtk_data;
  for (auto k = 0u; k < U_DIM; ++k) {
    vtk_data[k].type = T8_VTK_SCALAR;
    strcpy (vtk_data[k].description, ("u" + std::to_string (k)).c_str ());
    vtk_data[k].data = element_data[k];
  }

  const t8_element_t *element;
  const auto num_local_trees = t8_forest_get_num_local_trees (forest);

  auto current_index = 0u;
  for (auto tree_idx = 0u, current_index = 0u; tree_idx < num_local_trees; ++tree_idx) {
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    auto *data = t8_mra::get_mra_forest_data<T> (forest);
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
main (int argc, char **argv)
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
  // auto f3 = [] (double x, double y) -> std::array<double, 1> { return { 1.0 }; };

  auto max_level = 6u;
  auto c_thresh = 1.0;
  auto gamma = 1.0;  /// Order of convergence
  auto dunavant_rule = 10;

  bool balanced = false;

  constexpr int P = 3;
  constexpr int U = 1;

  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;
  using mra_type = t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P>;

  auto *test_scheme = t8_scheme_new_default ();
  t8_cmesh_t cmesh = t8_cmesh_square (comm);
  // t8_cmesh_t cmesh = t8_cmesh_l_shape (comm);
  printf ("created test_scheme and cmesh\n");

  mra_type mra_test (max_level, c_thresh, gamma, dunavant_rule, balanced, comm);
  printf ("created mra object\n");

  /// TODO f4 seems to be buggy
  mra_test.initialize_data (cmesh, test_scheme, max_level, f3);
  // printf ("create mra object\n");
  // mra_test.multiscale_transformation (0, max_level);
  // printf ("Size init data: %zu\n", mra_test.get_lmi_map ()->size ());
  // printf ("Initialize data\n");
  // mra_test.inverse_multiscale_transformation (0, max_level);
  // printf ("Size init data: %zu\n", mra_test.get_lmi_map ()->size ());

  printf ("Size init data: %zu\n", mra_test.get_lmi_map ()->size ());

  // Manual verification: Check first element
  {
    auto *mra_data = t8_mra::get_mra_forest_data<element_data_type> (mra_test.forest);
    const auto *element = t8_forest_get_leaf_element_in_tree (mra_test.forest, 0, 0);
    const auto lmi = t8_mra::get_lmi_from_forest_data<element_data_type> (mra_data, 0);
    const auto &elem_data = mra_data->lmi_map->get (lmi);
    const auto volume = t8_forest_element_volume (mra_test.forest, 0, element);

    printf ("\n=== Manual Verification (First Element) ===\n");
    printf ("Volume: %.8e\n", volume);
    printf ("Mean value: %.8e\n", t8_mra::mean_val<element_data_type> (mra_test.forest, 0, lmi, element)[0]);

    // Get physical vertices DIRECTLY from t8code (no permutation yet)
    printf ("T8code vertices (BEFORE permutation):\n");
    for (int v = 0; v < 3; ++v) {
      double coords[3];
      t8_forest_element_coordinate (mra_test.forest, 0, element, v, coords);
      printf ("  t8_vertex_%d: (%.6f, %.6f)\n", v, coords[0], coords[1]);
    }

    printf ("Vertex order from elem_data: [%d, %d, %d]\n", elem_data.order[0], elem_data.order[1], elem_data.order[2]);

    // Get physical vertices
    double vertices[9];
    for (int v = 0; v < 3; ++v) {
      double coords[3];
      t8_forest_element_coordinate (mra_test.forest, 0, element, v, coords);
      const int perm_v = elem_data.order[v];
      vertices[3 * perm_v] = coords[0];
      vertices[3 * perm_v + 1] = coords[1];
      vertices[3 * perm_v + 2] = coords[2];
    }

    printf ("Vertices (AFTER permutation for DG reference):\n");
    for (int v = 0; v < 3; ++v) {
      printf ("  v%d: (%.6f, %.6f)\n", v, vertices[3 * v], vertices[3 * v + 1]);
    }

    // TEST: Try swapping vertices 1 and 2 manually
    printf ("\n=== TESTING: Swap vertices 1 and 2 ===\n");
    double vertices_swapped[9];
    for (int v = 0; v < 3; ++v) {
      double coords[3];
      t8_forest_element_coordinate (mra_test.forest, 0, element, v, coords);
      // Map: t8_vertex_0 -> ref_vertex_0, t8_vertex_1 -> ref_vertex_2, t8_vertex_2 -> ref_vertex_1
      int ref_v = (v == 1) ? 2 : (v == 2) ? 1 : 0;
      vertices_swapped[3 * ref_v] = coords[0];
      vertices_swapped[3 * ref_v + 1] = coords[1];
      vertices_swapped[3 * ref_v + 2] = coords[2];
    }
    printf ("Swapped vertices:\n");
    for (int v = 0; v < 3; ++v) {
      printf ("  v%d: (%.6f, %.6f)\n", v, vertices_swapped[3 * v], vertices_swapped[3 * v + 1]);
    }

    // Evaluate at some test points
    printf ("DG evaluation at test points:\n");
    std::vector<std::array<double, 2>> test_pts = {
      { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.5, 0.0 }, { 0.5, 0.5 }, { 0.0, 0.5 }, { 1.0 / 3.0, 1.0 / 3.0 }
    };

    for (const auto &pt : test_pts) {
      auto val = t8_mra::evaluate_dg_at_point (elem_data, pt[0], pt[1], volume);
      // Map to physical coordinates
      const double lambda0 = 1.0 - pt[0] - pt[1];
      const double lambda1 = pt[0];
      const double lambda2 = pt[1];
      double phys_x = lambda0 * vertices[0] + lambda1 * vertices[3] + lambda2 * vertices[6];
      double phys_y = lambda0 * vertices[1] + lambda1 * vertices[4] + lambda2 * vertices[7];

      // Evaluate exact function at physical point
      auto exact = f3 (phys_x, phys_y);

      printf ("  ref(%.3f, %.3f) -> phys(%.6f, %.6f): DG=%.6e, exact=%.6e, diff=%.3e\n", pt[0], pt[1], phys_x, phys_y,
              val[0], exact[0], std::abs (val[0] - exact[0]));
    }
    printf ("==========================================\n\n");
  }
  // t8_write_vtu<element_data_type> (mra_test.forest,
  //                                  ("refinement_test/test_init_" + std::to_string (max_level)).c_str ());
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra_test.forest, "refinement_test/lagrange_init", P,
                                                        true);  // Quadratic Lagrange with debug

  mra_test.coarsening_new (0u, max_level);
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra_test.forest, "refinement_test/lagrange_coarse", P,
                                                        false);  // Quadratic Lagrange with debug
  // t8_write_vtu<element_data_type> (mra_test.forest,
  //                                  ("refinement_test/test_step_1_" + std::to_string (max_level)).c_str ());
  // printf ("Size first coarsening: %zu\n", mra_test.get_lmi_map ()->size ());
  //
  mra_test.refinement_new (0u, max_level);
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra_test.forest, "refinement_test/lagrange_refine", P,
                                                        false);  // Quadratic Lagrange with debug

  mra_test.coarsening_new (0u, max_level);
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra_test.forest, "refinement_test/lagrange_coarse2", P,
                                                        false);  // Quadratic Lagrange with debug
  // t8_write_vtu<element_data_type> (mra_test.forest,
  //                                  ("refinement_test/test_step_2_" + std::to_string (max_level)).c_str ());
  // printf ("Size first refinement: %zu\n", mra_test.get_lmi_map ()->size ());
  //
  // mra_test.coarsening_new (0u, max_level);
  // t8_write_vtu<element_data_type> (mra_test.forest,
  //                                  ("refinement_test/test_step_3_" + std::to_string (max_level)).c_str ());
  // printf ("Size second coarsening: %zu\n", mra_test.get_lmi_map ()->size ());

  // mra_test.coarsening (0, max_level);
  // // mra_test.refinement (1, 2);
  // printf ("Size after coarsening: %zu\n", mra_test.get_lmi_map ()->size ());
  // // mra_test.refinement (1, 2);
  // // printf ("Size after refinement: %zu\n", mra_test.get_lmi_map ()->size ());
  //
  // t8_write_vtu<element_data_type> (mra_test.forest,
  //                                  ("f3_refine_test" + std::to_string (P) + "_" + std::to_string (max_level)).c_str ());
  //
  // printf ("Finished writing file\n");

  mra_test.cleanup ();
  printf ("cleaned everything...\n");

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  printf ("Program finished...\n");
  return 0;
}
