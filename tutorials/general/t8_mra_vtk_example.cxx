/**
 * Example demonstrating VTK output for MRA with higher-order Lagrange cells
 *
 * This example shows how to:
 * 1. Initialize a uniform forest with modal DG data
 * 2. Export the data using higher-order Lagrange cells for VTK
 * 3. Compare cell-averaged output vs. full polynomial representation
 */

#include <t8_mra/t8_mra.hpp>
#include <t8_mra/t8_mra_vtk.hpp>
#include <t8_mra/t8_mra_vtk_subdivide.hpp>

#include "t8.h"
#include "t8_cmesh.hxx"
#include "t8_eclass.h"
#include "t8_forest/t8_forest_general.h"
#include "t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx"
#include "t8_mra/data/cell_data.hpp"

// Simple test function: quarter circle with jump
auto test_function = [] (double x, double y) -> std::array<double, 1> {
  double r = x * x + y * y;
  return { (r < 0.25) ? (x * y + x + 3.) : (x * x * y - 2. * x * y * y + 3. * x) };
};

// Smooth test function: Gaussian blob
auto gaussian_blob = [] (double x, double y) -> std::array<double, 1> {
  double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
  return { std::exp (-50.0 * r2) };
};

// Create a simple square mesh
t8_cmesh_t
create_square_mesh (sc_MPI_Comm comm)
{
  double vertices[18] = {
    0, 0, 0, 0, 1, 0, 1, 0, 0,  //triangle 1
    1, 1, 0, 1, 0, 0, 0, 1, 0,  //triangle 2
  };

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);

  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);

  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);
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

  // MRA parameters
  auto max_level = 6u;
  auto c_thresh = 1.0;
  auto gamma = 1.0;
  auto dunavant_rule = 10;
  bool balanced = false;

  // Polynomial order for DG basis
  constexpr int P = 3;  // Cubic polynomials
  constexpr int U = 1;  // Scalar function

  using element_data_type = t8_mra::data_per_element<T8_ECLASS_TRIANGLE, U, P>;
  using mra_type = t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P>;

  printf ("=== MRA VTK Output Example ===\n");
  printf ("Polynomial order P = %d\n", P);
  printf ("Number of DOF per element = %d\n", (int) element_data_type::DOF);

  // Create mesh and MRA object
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_t cmesh = create_square_mesh (comm);

  mra_type mra (max_level, c_thresh, gamma, dunavant_rule, balanced, comm);

  // Initialize with a test function at a specific level
  int init_level = 6;
  printf ("\nInitializing uniform forest at level %d...\n", init_level);
  mra.initialize_data (cmesh, scheme, init_level, test_function);

  printf ("Total elements: %zu\n", mra.get_lmi_map ()->size ());

  // ===== METHOD 1: Cell-averaged values (fast, simple) =====
  printf ("\n--- Writing cell-averaged VTK output ---\n");
  printf ("This shows only the mean value per cell (P0 representation)\n");
  t8_mra::write_forest_cell_average_vtk<element_data_type> (mra.forest, "mra_output_cell_average");

  // ===== METHOD 2: Higher-order Lagrange representation =====
  printf ("\n--- Writing Lagrange P2 (quadratic) VTK output (with debug) ---\n");
  printf ("This evaluates the DG polynomial at Lagrange nodes for higher-order visualization\n");
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.forest, "mra_output_lagrange_p2", 2,
                                                        true);  // Quadratic Lagrange with debug

  printf ("\n--- Writing Lagrange P3 (cubic) VTK output ---\n");
  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.forest, "mra_output_lagrange_p3", 3,
                                                        false);  // Cubic Lagrange

  // ===== METHOD 3: Subdivided mesh (guaranteed to work in ParaView) =====
  printf ("\n--- Writing subdivided VTK output (16x16 subdivisions) ---\n");
  printf ("This subdivides each element into linear triangles for better ParaView compatibility\n");
  printf ("For polynomial order P=%d, using %d subdivisions per edge\n", P, 16);
  t8_mra::write_forest_subdivided_vtk<element_data_type> (mra.forest, "mra_output_subdivided", 16);

  // ===== Optional: Test with adapted mesh =====
  printf ("\n--- Testing with adaptive refinement ---\n");
  mra.coarsening_new (0u, max_level);
  printf ("After coarsening: %zu elements\n", mra.get_lmi_map ()->size ());

  t8_mra::write_forest_lagrange_vtk<element_data_type> (mra.forest, "mra_output_adapted", 2);

  // Cleanup
  mra.cleanup ();

  printf ("\n=== VTK files written successfully ===\n");
  printf ("Open with ParaView or VisIt for visualization\n");
  printf ("Compare 'cell_average' (P0) with 'lagrange_p2/p3' (higher-order) output\n");

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
