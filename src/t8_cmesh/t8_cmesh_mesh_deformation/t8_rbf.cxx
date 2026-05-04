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

typedef Eigen::SparseMatrix<double> SpMat;  // declares a column-major sparse matrix type of double

/**
   * Solves the linear system A * alpha = displacements to find the weight. 
   * This step is mandatory to later be able to interpolate the inner nodes.
   */
void
t8_rbf::solve ()
{
  double estimation_of_entries
    = 50;  // brauchen noch gute Schätzung für die Anzahl der nicht-Null Element in der Matrix A für compactly supported RBFs
  const size_t num_boundary_nodes = boundary_nodes.size ();
  /** Check if there are any boundary nodes. */
  if (num_boundary_nodes == 0) {
    t8_errorf ("ERROR: Boundary nodes are not added correctly.\n.");
    SC_ABORTF ("The current RBF instance has no boundary nodes.");
  }

  /** Fill the displacement vector with the values of the boundary nodes' displacements. */
  Eigen::MatrixXd displacements (num_boundary_nodes, 3);
  for (size_t i = 0; i < num_boundary_nodes; ++i) {
    displacements.row (i) = Eigen::Vector3d (boundary_nodes[i].displacement[0], boundary_nodes[i].displacement[1],
                                             boundary_nodes[i].displacement[2]);
  }
  /** Allocate the memory for the weight Vector.  */
  Eigen::MatrixXd alpha (num_boundary_nodes, 3);

  if (rbf_function->is_compactly_supported ()) {
    /** If the RBF used is compactly supported use a sparse matrix. Then the conjugate gradient method can be used for solving the linear equation system. */
    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> coefficients;
    coefficients.reserve (num_boundary_nodes * estimation_of_entries);
    /** Fill the matrix A with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
    for (size_t row = 0; row < num_boundary_nodes; ++row) {
      /** Because of the symmetric property of the distance between nodes, we only need to compute the upper triangular part of the matrix 
       * and can mirror the values to the lower triangular part. */
      for (size_t col = row; col < num_boundary_nodes; ++col) {
        /** Calculate the distance between the current pair of boundary nodes. */
        const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
        /** Evaluate the radial basis function for the current distance. */
        const double psi = rbf_function->evaluate (distance);
        coefficients.emplace_back (row, col, psi);
        if (row != col) {
          coefficients.emplace_back (col, row, psi);
        }
      }
    }
    Eigen::SparseMatrix<double> A (num_boundary_nodes, num_boundary_nodes);
    A.setFromTriplets (coefficients.begin (), coefficients.end ());
    Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> solver;
    solver.compute (A);
    if (solver.info () != Eigen::Success) {
      t8_errorf ("ERROR: Decomposition of the matrix A failed.\n.");
      SC_ABORTF ("The linear system cannot be solved.");
    }
  }
  else {
    /** If the RBF used is globally supported use a dense matrix. */
  }

  /** Create the matrix A for the linear system. */
  std::vector<double> A (num_boundary_nodes * num_boundary_nodes);
  /** Fill the matrix A with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
  for (size_t row = 0; row < num_boundary_nodes; ++row) {
    /** Because of the symmetric property of the distance between nodes, we only need to compute the upper triangular part of the matrix 
       * and can mirror the values to the lower triangular part. */
    for (size_t col = row; col < num_boundary_nodes; ++col) {
      /** Calculate the distance between the current pair of boundary nodes. */
      const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
      /** Evaluate the radial basis function for the current distance. */
      const double psi = rbf_function->evaluate (distance);
      /** Write the value to the matrix A. */
      A[row * num_boundary_nodes + col] = psi;
      /** Mirror the value to the lower triangular part of the matrix if it's not on the diagonal. */
      if (col != row) {
        A[col * num_boundary_nodes + row] = psi;
      }
    }
  }
#if 0
  /** Create the right-hand side vector for the linear system. It contains the displacements of the boundary nodes. */
  std::vector<double> displacements (num_boundary_nodes * 3);
  for (size_t i = 0; i < num_boundary_nodes; ++i) {
    displacements[3 * i + 0] = boundary_nodes[i].displacement[0];
    displacements[3 * i + 1] = boundary_nodes[i].displacement[1];
    displacements[3 * i + 2] = boundary_nodes[i].displacement[2];
  }
#endif
}
