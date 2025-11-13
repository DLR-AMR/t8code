#include <t8_mra/t8_mra.hpp>
#include <t8_mra/t8_mra_vtk.hpp>

#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_mra/data/cell_data.hpp"

// Simple test function: constant function f(x,y) = 1.0
auto constant_function = [] (double x, double y) -> std::array<double, 1> {
  return { 1.0 };
};

// Linear function f(x,y) = x + y
auto linear_function = [] (double x, double y) -> std::array<double, 1> {
  return { x + y };
};

t8_cmesh_t
create_single_triangle (sc_MPI_Comm comm)
{
  double vertices[9] = {
    0, 0, 0,  // vertex 0
    1, 0, 0,  // vertex 1
    0, 1, 0,  // vertex 2
  };

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
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

  constexpr int P = 3;
  constexpr int U = 1;
  auto max_level = 5u;

  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;
  using mra_type = t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P>;

  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_t cmesh = create_single_triangle (comm);

  mra_type mra (max_level, 1.0, 1.0, 10, false, comm);

  printf ("=== Testing with constant function f(x,y) = 1.0 ===\n");
  mra.initialize_data (cmesh, scheme, 0, constant_function);

  // Get the data
  auto *mra_data = t8_mra::get_mra_forest_data<element_data_type> (mra.forest);
  auto lmi = t8_mra::get_lmi_from_forest_data<element_data_type> (mra_data, 0);
  const auto &elem_data = mra_data->lmi_map->get (lmi);

  printf ("Number of DOF: %d\n", (int)element_data_type::DOF);
  printf ("Stored coefficients:\n");
  for (auto i = 0u; i < element_data_type::DOF; ++i) {
    printf ("  u_coeffs[%u] = %.10e\n", i, elem_data.u_coeffs[element_data_type::dg_idx (0, i)]);
  }

  // Get volume
  const auto *element = t8_forest_get_leaf_element_in_tree (mra.forest, 0, 0);
  const auto volume = t8_forest_element_volume (mra.forest, 0, element);
  printf ("Element volume: %.10e\n", volume);

  // Test mean value function
  const auto mean = t8_mra::mean_val<element_data_type> (mra.forest, 0, lmi, element);
  printf ("Mean value (from mean_val function): %.10e\n", mean[0]);

  // Test basis evaluation at different points
  printf ("\n=== Testing basis function evaluation ===\n");
  for (auto i = 0u; i < element_data_type::DOF; ++i) {
    double phi_at_origin = t8_mra::skalierungsfunktion (i, 0.0, 0.0);
    double phi_at_center = t8_mra::skalierungsfunktion (i, 1.0/3.0, 1.0/3.0);
    printf ("  phi_%u(0,0) = %.10e, phi_%u(1/3,1/3) = %.10e\n",
            i, phi_at_origin, i, phi_at_center);
  }

  // Test evaluation at several points
  printf ("\n=== Testing function evaluation at points ===\n");
  std::vector<std::array<double, 2>> test_points = {
    {0.0, 0.0},      // corner
    {1.0, 0.0},      // corner
    {0.0, 1.0},      // corner
    {1.0/3.0, 1.0/3.0},  // center
    {0.5, 0.0},      // edge midpoint
  };

  for (const auto &pt : test_points) {
    auto eval = t8_mra::evaluate_dg_at_point (elem_data, pt[0], pt[1], volume);
    printf ("  f(%.3f, %.3f) = %.10e (expected: 1.0)\n", pt[0], pt[1], eval[0]);
  }

  mra.cleanup ();

  printf ("\n=== Testing with linear function f(x,y) = x + y ===\n");
  t8_cmesh_t cmesh2 = create_single_triangle (comm);
  mra_type mra2 (max_level, 1.0, 1.0, 10, false, comm);
  mra2.initialize_data (cmesh2, scheme, 0, linear_function);

  auto *mra_data2 = t8_mra::get_mra_forest_data<element_data_type> (mra2.forest);
  auto lmi2 = t8_mra::get_lmi_from_forest_data<element_data_type> (mra_data2, 0);
  const auto &elem_data2 = mra_data2->lmi_map->get (lmi2);

  printf ("Stored coefficients:\n");
  for (auto i = 0u; i < element_data_type::DOF; ++i) {
    printf ("  u_coeffs[%u] = %.10e\n", i, elem_data2.u_coeffs[element_data_type::dg_idx (0, i)]);
  }

  const auto *element2 = t8_forest_get_leaf_element_in_tree (mra2.forest, 0, 0);
  const auto volume2 = t8_forest_element_volume (mra2.forest, 0, element2);

  printf ("\nEvaluating at points:\n");
  for (const auto &pt : test_points) {
    auto eval = t8_mra::evaluate_dg_at_point (elem_data2, pt[0], pt[1], volume2);
    double expected = pt[0] + pt[1];
    printf ("  f(%.3f, %.3f) = %.10e (expected: %.3f)\n", pt[0], pt[1], eval[0], expected);
  }

  mra2.cleanup ();

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
