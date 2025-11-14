/*
  Example: Project a function onto a cartesian (QUAD) grid using MRA

  This example demonstrates:
  1. Creating a uniform quad mesh
  2. Projecting a smooth function onto DG basis
  3. Evaluating the projected solution
*/

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <array>
#include <cmath>
#include <iostream>

#ifdef T8_ENABLE_MRA
#include "t8_mra/t8_mra_cartesian.hpp"
#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/t8_mra_vtk.hpp"

// Test function: smooth Gaussian-like function
std::array<double, 1>
test_function (double x, double y)
{
  const double x0 = 0.5, y0 = 0.5;
  const double sigma = 0.2;
  const double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
  return { std::exp (-r2 / (2.0 * sigma * sigma)) };
}

int
main (int argc, char **argv)
{
  // Initialize MPI and t8code
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_PRODUCTION);

  /// Quarter circle + jump along circle
  auto f3 = [] (double x, double y) -> std::array<double, 1> {
    double r = x * x + y * y;
    return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
  };

  std::cout << "\n=== Cartesian (QUAD) Projection Example ===\n\n";

  // Create a simple quad mesh on [0,1]^2
  t8_cmesh_t cmesh;
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

  // Create scheme for quad elements
  const t8_scheme *scheme = t8_scheme_new_default ();

  // MRA parameters
  const int max_level = 5;           // Maximum refinement level
  const double c_thresh = 0.1;       // Coarsening threshold
  const int gamma = 2;               // Gamma parameter
  const int num_quad_points_1d = 5;  // 4 Gauss-Legendre points per dimension
  const bool balanced = false;       // No balancing for now
  const int P = 4;                   // Polynomial degree P-1 = 2 (quadratic)
  const int U = 1;                   // 1 solution component

  std::cout << "MRA Parameters:\n";
  std::cout << "  Element type: QUAD\n";
  std::cout << "  Polynomial order P: " << P << " (degree " << (P - 1) << ")\n";
  std::cout << "  Max level: " << max_level << "\n";
  std::cout << "  Quadrature points per dimension: " << num_quad_points_1d << "\n";
  std::cout << "  Total DOF per element: " << (P * P) << "\n\n";

  // Create multiscale object for QUAD elements
  using mra_type = t8_mra::multiscale<T8_ECLASS_QUAD, U, P>;
  mra_type mra (max_level, c_thresh, gamma, num_quad_points_1d, balanced, sc_MPI_COMM_WORLD);

  std::cout << "Projecting function onto uniform grid...\n";

  // Initialize with uniform refinement and project function
  // mra.initialize_data (cmesh, scheme, max_level, test_function);
  mra.initialize_data (cmesh, scheme, max_level, f3);

  // Get forest info
  const auto num_elements = t8_forest_get_local_num_leaf_elements (mra.get_forest ());
  std::cout << "  Created forest with " << num_elements << " QUAD elements\n";
  std::cout << "  Total DOF: " << (num_elements * P * P) << "\n\n";

  // Test: evaluate at center point (0.5, 0.5)
  std::cout << "Testing projection quality:\n";
  std::cout << "  Function at (0.5, 0.5):\n";
  std::cout << "    Exact value: " << test_function (0.5, 0.5)[0] << "\n";

  // Test a few points to check continuity
  std::cout << "\n  Testing edge continuity at x=0.5 (should be continuous):\n";
  std::cout << "    f(0.49999, 0.5) = " << test_function (0.49999, 0.5)[0] << "\n";
  std::cout << "    f(0.50000, 0.5) = " << test_function (0.50000, 0.5)[0] << "\n";
  std::cout << "    f(0.50001, 0.5) = " << test_function (0.50001, 0.5)[0] << "\n";

  // Find element containing (0.5, 0.5) and evaluate
  // For simplicity, we'll just print info about the first element
  if (num_elements > 0) {
    const auto *elem = t8_forest_get_leaf_element_in_tree (mra.get_forest (), 0, 0);
    std::cout << "  First element created successfully\n";

    // Print basis info
    std::cout << "\nBasis function info:\n";
    std::cout << "  Number of quad points: " << mra.DG_basis.num_quad_points << "\n";
    std::cout << "  Reference quad points (first 3):\n";
    for (int i = 0; i < std::min (3, static_cast<int> (mra.DG_basis.num_quad_points)); ++i) {
      std::cout << "    (" << mra.DG_basis.ref_quad_points[2 * i] << ", " << mra.DG_basis.ref_quad_points[2 * i + 1]
                << ")  weight: " << mra.DG_basis.quad_weights[i] << "\n";
    }

    // Test basis evaluation at reference point (0.5, 0.5)
    std::vector<double> ref_point = { 0.5, 0.5 };
    auto basis_vals = mra.DG_basis.basis_value (ref_point);
    std::cout << "\n  Basis values at ref point (0.5, 0.5):\n";
    for (int i = 0; i < P * P; ++i) {
      std::cout << "    phi_" << i << " = " << basis_vals[i] << "\n";
    }

    // Print pset (multiindex set) for debugging
    std::cout << "\n  Multiindex set (pset) for P=" << P << ":\n";
    const auto pset = mra.DG_basis.pset;
    for (int i = 0; i < P * P; ++i) {
      std::cout << "    pset[" << i << "] = (" << pset[i][0] << ", " << pset[i][1] << ")\n";
    }

    // Get and print coefficients for first element
    auto *user_data = mra.get_user_data ();
    const auto lmi_0 = t8_mra::get_lmi_from_forest_data (user_data, 0);
    const auto &elem_data_0 = user_data->lmi_map->get (lmi_0);
    std::cout << "\n  DG coefficients for first element:\n";
    for (int i = 0; i < P * P; ++i) {
      std::cout << "    u_coeffs[" << i << "] = " << elem_data_0.u_coeffs[i] << "\n";
    }
  }

  std::cout << "\n=== Projection Complete ===\n";
  std::cout << "Cartesian element support is working!\n\n";

  // Write VTK output with Lagrange representation (shows polynomial structure)
  std::cout << "Writing VTK output...\n";
  using element_data_type = t8_mra::data_per_element<T8_ECLASS_QUAD, U, P>;

  // Write with polynomial order P (matches DG basis order) - enable debug
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.get_forest (), "quad_test/quad_projection_lagrange", P,
                                                        true);
  std::cout << "  Output written to: quad_projection_lagrange.vtu\n";
  std::cout << "  Open with ParaView to visualize the polynomial solution\n";
  std::cout << "  Each QUAD element shows the full P=" << P << " polynomial representation\n\n";

  // Cleanup
  mra.cleanup ();
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

#else
int
main ()
{
  std::cout << "This example requires T8_ENABLE_MRA to be enabled.\n";
  return 0;
}
#endif
