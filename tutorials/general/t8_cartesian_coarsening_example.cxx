/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_cartesian_coarsening_example.cxx
 * Demonstrates adaptive coarsening for cartesian (QUAD) elements using multiscale decomposition.
 *
 * This example:
 * 1. Creates a uniform QUAD mesh at a fine level
 * 2. Projects a smooth Gaussian function onto it
 * 3. Performs wavelet-based adaptive coarsening
 * 4. Writes VTK output showing the adapted mesh
 *
 * The coarsening uses detail coefficients from multiscale transformation to identify
 * regions where the solution is smooth (small details) and can be coarsened.
 */

#include <iostream>
#include <cmath>
#include <array>

#include <t8.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>

#include <t8_mra/t8_mra_cartesian.hpp>
#include <t8_mra/t8_mra_vtk.hpp>

/**
 * @brief Gaussian blob function centered at (0.5, 0.5)
 *
 * @param x X-coordinate
 * @param y Y-coordinate
 * @return std::array<double, 1> Function value (single component)
 */
std::array<double, 1>
gaussian_function (double x, double y)
{
  const double x0 = 0.5;     // Center x
  const double y0 = 0.5;     // Center y
  const double sigma = 0.2;  // Standard deviation

  const double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
  const double val = std::exp (-r2 / (2.0 * sigma * sigma));

  return { val };
  // return { 1.0 };
}

int
main (int argc, char **argv)
{
  // Initialize MPI and t8code
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, nullptr, SC_LP_ESSENTIAL);
  t8_init (SC_LP_ESSENTIAL);

  std::cout << "=== Cartesian (QUAD) Adaptive Coarsening Example ===\n\n";

  // Test polynomial f(x,y) = x*y (degree 2, should be exactly representable with P=3)
  // auto f3 = [] (double x, double y) -> std::array<double, 1> { return { 1.0 }; };
  // auto f3 = [] (double x, double y) -> std::array<double, 1> { return { x * y }; };
  // auto f3 = [] (double x, double y) -> std::array<double, 1> { return { x }; };
  auto f3 = [] (double x, double y) -> std::array<double, 1> {
    double r = x * x + y * y;
    return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
  };

  // MRA parameters
  constexpr int U = 1;                  // Number of solution components
  constexpr int P = 3;                  // Polynomial order (degree = P-1 = 2)
  const int max_level = 6;              // Maximum level for adaptation
  const int min_level = 0;              // Allow coarsening down to level 2
  const int initial_level = max_level;  // Start with 2^6 x 2^6 = 4096 elements

  // Coarsening parameters
  const double c_thresh = 1.0;       // Threshold for coarsening (smaller = more aggressive)
  const int gamma = 1;               // Parameter for threshold scaling
  const int num_quad_points_1d = 4;  // Gauss-Legendre quadrature points per dimension
  const bool balanced = false;       // Whether to balance the forest

  using element_data_type = t8_mra::data_per_element<T8_ECLASS_QUAD, U, P>;

  std::cout << "MRA Parameters:\n";
  std::cout << "  Element type: QUAD\n";
  std::cout << "  Polynomial order P: " << P << " (degree " << (P - 1) << ")\n";
  std::cout << "  Initial level: " << initial_level << " (" << (1 << initial_level) << "x" << (1 << initial_level)
            << " = " << ((1 << initial_level) * (1 << initial_level)) << " elements)\n";
  std::cout << "  Coarsening threshold: " << c_thresh << "\n";
  std::cout << "  Gamma parameter: " << gamma << "\n";
  std::cout << "  Quadrature points per dimension: " << num_quad_points_1d << "\n";
  std::cout << "  Total DOF per element: " << element_data_type::DOF << "\n\n";

  // Create multiscale object
  t8_mra::multiscale<T8_ECLASS_QUAD, U, P> mra (max_level, c_thresh, gamma, num_quad_points_1d, balanced,
                                                sc_MPI_COMM_WORLD);

  // Create unit square mesh
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();

  std::cout << "Projecting Gaussian function onto uniform level " << initial_level << " grid...\n";

  // Initialize with uniform mesh and project function
  // mra.initialize_data (cmesh, scheme, initial_level, gaussian_function);
  mra.initialize_data (cmesh, scheme, initial_level, f3);

  const auto num_elements_initial = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  std::cout << "  Created forest with " << num_elements_initial << " QUAD elements\n";
  std::cout << "  Total DOF: " << (num_elements_initial * element_data_type::DOF) << "\n\n";

  // Write initial uniform mesh
  std::cout << "Writing initial uniform mesh VTK output...\n";
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.get_forest (), "quad_coarsening/initial_uniform", P - 1,
                                                        true);
  std::cout << "  Output written to: quad_coarsening/initial_uniform.vtu\n\n";

  // Perform adaptive coarsening
  std::cout << "=== Starting Adaptive Coarsening ===\n";
  std::cout << "Coarsening from level " << max_level << " down to level " << min_level << "...\n\n";

  mra.coarsening_new (min_level, max_level);

  const auto num_elements_final = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
  const double reduction_percent = 100.0 * (1.0 - static_cast<double> (num_elements_final) / num_elements_initial);

  std::cout << "\n=== Coarsening Complete ===\n";
  std::cout << "Initial elements: " << num_elements_initial << "\n";
  std::cout << "Final elements:   " << num_elements_final << "\n";
  std::cout << "Reduction:        " << reduction_percent << "%\n";
  std::cout << "Final DOF:        " << (num_elements_final * element_data_type::DOF) << "\n\n";

  // Write adapted mesh
  std::cout << "Writing adapted mesh VTK output...\n";
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.get_forest (), "quad_coarsening/adapted_mesh", P - 1, true);
  // t8_mra::write_forest_cell_average_vtk<element_data_type> (mra.get_forest (), "quad_coarsening/adapted_mesh");
  std::cout << "  Output written to: quad_coarsening/adapted_mesh.vtu\n\n";

  std::cout << "=== Visualization Instructions ===\n";
  std::cout << "Open the VTK files in ParaView:\n";
  std::cout << "  1. Load both 'initial_uniform.vtu' and 'adapted_mesh.vtu'\n";
  std::cout << "  2. Apply 'Surface' representation\n";
  std::cout << "  3. Color by 'u0' to see the solution\n";
  std::cout << "  4. Compare element sizes - adapted mesh has larger elements in smooth regions\n";
  std::cout << "  5. Enable 'Surface With Edges' to see the adaptive grid\n\n";

  std::cout << "Expected behavior:\n";
  std::cout << "  - Fine elements near the Gaussian peak (high gradient)\n";
  std::cout << "  - Coarse elements far from center (smooth, nearly zero)\n";
  std::cout << "  - Wavelet-based coarsening preserves solution accuracy\n\n";

  // Cleanup
  mra.cleanup ();
  // t8_scheme_unref (&scheme);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  std::cout << "Example completed successfully!\n";

  return 0;
}
