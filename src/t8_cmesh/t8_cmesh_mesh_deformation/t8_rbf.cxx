/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_mesh_deformation_rbf.cxx
 *  This file implements the Radial Basis Functions for the mesh deformation.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_rbf.hxx>
#include <t8_types/t8_vec.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

/**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
void
t8_rbf::solve ()
{
  const size_t num_boundary_nodes = boundary_nodes.size ();
  /** Check if there are any boundary nodes. */
  if (num_boundary_nodes == 0) {
    t8_errorf ("ERROR: Boundary nodes are not added correctly. The current RBF instance has no boundary nodes\n.");
  }

  /** Fill the displacement vector with the values of the boundary nodes' displacements. */
  Eigen::MatrixXd displacements (num_boundary_nodes, 3);
  for (size_t i = 0; i < num_boundary_nodes; ++i) {
    displacements.row (i) = Eigen::Vector3d (boundary_nodes[i].displacement[0], boundary_nodes[i].displacement[1],
                                             boundary_nodes[i].displacement[2]);
  }
  /** Weight vector.  */
  Eigen::MatrixXd alpha;

  if (rbf_function->is_compactly_supported ()) {
    alpha = solve_compactly_supported_rbf (displacements, num_boundary_nodes);
  }
  else {
    alpha = solve_globally_supported_rbf (displacements, num_boundary_nodes);
  }
  /** Copy the calculated weights back to the boundary nodes. */
  for (size_t i = 0; i < num_boundary_nodes; ++i) {
    boundary_nodes[i].weight[0] = alpha (i, 0);
    boundary_nodes[i].weight[1] = alpha (i, 1);
    boundary_nodes[i].weight[2] = alpha (i, 2);
  }
}

Eigen::MatrixXd
t8_rbf::solve_compactly_supported_rbf (const Eigen::MatrixXd &displacements, const size_t num_boundary_nodes) const
{
  /** The RBF used is compactly supported so we can use a sparse matrix. The conjugate gradient method can be used for solving the linear equation system. */
  size_t estimation_of_entries
    = 50;  // brauchen noch gute Schätzung für die Anzahl der nicht-Null Element in der Matrix A für compactly supported RBFs
  typedef Eigen::Triplet<double> triplet;
  std::vector<triplet> coefficients;
  coefficients.reserve (num_boundary_nodes * estimation_of_entries);
  /** Fill the matrix A with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
  for (size_t row = 0; row < num_boundary_nodes; ++row) {

    for (size_t col = 0; col < num_boundary_nodes; ++col) {
      /** Calculate the distance between the current pair of boundary nodes. */
      const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
      /** Evaluate the radial basis function for the current distance. */
      const double psi = rbf_function->evaluate (distance, boundary_nodes[col].local_support_radius);
      /** Check if the basis function value (the influence factor) is not equal to zero with a numerical tolerance of 1e-12. */
      if (std::abs (psi) > 1e-12) {
        coefficients.emplace_back (row, col, psi);
      }
    }
  }
  /** Fill the sparse matrix A with the triplets. */
  Eigen::SparseMatrix<double> A (num_boundary_nodes, num_boundary_nodes);
  A.setFromTriplets (coefficients.begin (), coefficients.end ());
  /** Solve the linear system using the conjugate gradient method. */
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.setTolerance (1e-10);
  solver.compute (A);

  if (solver.info () != Eigen::Success) {
    t8_errorf ("ERROR: Decomposition of the matrix A failed. The linear system cannot be solved.\n");
  }

  return solver.solve (displacements);
  ;
}

Eigen::MatrixXd
t8_rbf::solve_globally_supported_rbf (const Eigen::MatrixXd &displacements, const size_t num_boundary_nodes) const
{
  /** The RBF used is globally supported so we need to use a dense matrix. The linear equation system can be solved with the built-in solver of Eigen. */

  /** Create the matrix A for the linear system. */
  Eigen::MatrixXd A (num_boundary_nodes, num_boundary_nodes);
  /** Fill the matrix A with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
  for (size_t row = 0; row < num_boundary_nodes; ++row) {
    /** Because of the symmetric property of the distance between nodes, we only need to compute the upper triangular part of the matrix 
       * and can mirror the values to the lower triangular part. */
    for (size_t col = row; col < num_boundary_nodes; ++col) {
      /** Calculate the distance between the current pair of boundary nodes. */
      const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
      /** Evaluate the radial basis function for the current distance. */
      const double psi = rbf_function->evaluate (distance, 0.0);
      /** Write the value to the matrix A. */
      A (row, col) = psi;
      /** Mirror the value to the lower triangular part of the matrix if it's not on the diagonal. */
      if (col != row) {
        A (col, row) = psi;
      }
    }
  }
  Eigen::LDLT<Eigen::MatrixXd> solver (A);

  if (solver.info () != Eigen::Success) {
    t8_errorf ("ERROR: The global RBF system could not be solved with LDLT.\n");
  }

  return solver.solve (displacements);
}
