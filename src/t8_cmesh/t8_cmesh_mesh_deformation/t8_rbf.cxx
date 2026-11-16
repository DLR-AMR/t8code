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

/** \file t8_rbf.cxx
 *  This file implements the Radial Basis Functions for the mesh deformation.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh_mesh_deformation/t8_rbf.hxx>
#include <t8_vector/t8_vec.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <limits>

/** Solve the RBF system for a compactly supported RBF.
 * The interpolation matrix is sparse, so a sparse direct solver is used.
 * \param [in] boundary_nodes The boundary nodes of the RBF system.
 * \param [in] rbf_function The compactly supported radial basis function.
 * \param [in] displacements The matrix of the boundary node displacements. Each row corresponds to a
 *                          boundary node and the three columns correspond to the x, y, and z components of the displacement.
 * \return The matrix of the calculated weights alpha.
 */
static Eigen::MatrixXd
t8_rbf_solve_compactly_supported (const std::vector<t8_rbf_boundary_node> &boundary_nodes,
                                  const t8_rbf_function &rbf_function, const Eigen::MatrixXd &displacements)
{
  const size_t num_boundary_nodes = boundary_nodes.size ();
  /* The RBF used is compactly supported so we can use a sparse matrix and a sparse direct solver. */
  typedef Eigen::Triplet<double> triplet;
  std::vector<triplet> coefficients;
  /* Each row has at least one non-zero element on the diagonal.
   * This only reserves memory for the triplet list. More entries are added automatically. */
  coefficients.reserve (num_boundary_nodes);
  /* Fill the interpolation matrix with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
  for (size_t row = 0; row < num_boundary_nodes; ++row) {

    for (size_t col = 0; col < num_boundary_nodes; ++col) {
      /* Calculate the distance between the current pair of boundary nodes. */
      const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
      /* Evaluate the radial basis function for the current distance. */
      const double psi = rbf_function.evaluate (distance, boundary_nodes[col].local_support_radius);
      /* Check if the basis function value (the influence factor) is not equal to zero with a numerical tolerance of 1e-12. */
      if (std::abs (psi) > 1e-12) {
        coefficients.emplace_back (row, col, psi);
      }
    }
  }
  /* Fill the sparse interpolation matrix with the triplets. */
  Eigen::SparseMatrix<double> interpolation_matrix (num_boundary_nodes, num_boundary_nodes);
  interpolation_matrix.setFromTriplets (coefficients.begin (), coefficients.end ());
  interpolation_matrix.makeCompressed ();
  /* Solve the linear system using the sparse LU decomposition. */
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute (interpolation_matrix);

  if (solver.info () != Eigen::Success) {
    t8_errorf ("ERROR: Decomposition of the interpolation matrix failed. The linear system cannot be solved.\n");
  }

  return solver.solve (displacements);
}

/** Solve the RBF system for a globally supported RBF.
 * The interpolation matrix is dense, so a dense solver is used.
 * \param [in] boundary_nodes The boundary nodes of the RBF system.
 * \param [in] rbf_function The globally supported radial basis function.
 * \param [in] displacements The matrix of the boundary node displacements. Each row corresponds to a
 *                          boundary node and the three columns correspond to the x, y, and z components of the displacement.
 * \return The matrix of the calculated weights alpha.
 */
static Eigen::MatrixXd
t8_rbf_solve_globally_supported (const std::vector<t8_rbf_boundary_node> &boundary_nodes,
                                 const t8_rbf_function &rbf_function, const Eigen::MatrixXd &displacements)
{
  const size_t num_boundary_nodes = boundary_nodes.size ();
  /* The RBF used is globally supported so we need to use a dense matrix. The linear equation system can be solved with the built-in solver of Eigen. */

  /* Create the interpolation matrix for the linear system. */
  Eigen::MatrixXd interpolation_matrix (num_boundary_nodes, num_boundary_nodes);
  /* Fill the interpolation matrix with the values of the radial basis function psi evaluated at the pairwise euclidean distances between all boundary nodes. */
  for (size_t row = 0; row < num_boundary_nodes; ++row) {
    /* Because of the symmetric property of the distance between nodes, we only need to compute the upper triangular part (including the diagonal) of the matrix
     * and can mirror the values to the lower triangular part. */
    for (size_t col = row; col < num_boundary_nodes; ++col) {
      /* Calculate the distance between the current pair of boundary nodes. */
      const double distance = t8_dist (boundary_nodes[row].position, boundary_nodes[col].position);
      /* Evaluate the radial basis function for the current distance and write the value to the interpolation matrix. */
      const double psi = rbf_function.evaluate (distance);
      interpolation_matrix (row, col) = psi;
      interpolation_matrix (col, row) = psi;
    }
  }
  Eigen::PartialPivLU<Eigen::MatrixXd> solver (interpolation_matrix);

  /* PartialPivLU assumes an invertible matrix, so we check the estimated reciprocal condition number. */
  if (solver.rcond () < std::numeric_limits<double>::epsilon ()) {
    t8_errorf ("ERROR: The interpolation matrix is singular or nearly singular. The linear system cannot be solved.\n");
  }

  return solver.solve (displacements);
}

/*
 * Solves the linear system A * alpha = displacements to find the weight.
 * This step is mandatory to later be able to interpolate the inner nodes.
 */
void
t8_rbf::solve ()
{
  const size_t num_boundary_nodes = boundary_nodes.size ();
  /* Check if there are any boundary nodes. */
  if (num_boundary_nodes == 0) {
    t8_errorf ("ERROR: Boundary nodes are not added correctly. The current RBF instance has no boundary nodes.\n");
    SC_ABORT ("The RBF instance has no boundary nodes.");
  }

  /* Fill the displacement vector with the values of the boundary nodes' displacements. */
  Eigen::MatrixXd displacements (num_boundary_nodes, 3);
  for (size_t inode = 0; inode < num_boundary_nodes; ++inode) {
    displacements.row (inode)
      = Eigen::Vector3d (boundary_nodes[inode].displacement[0], boundary_nodes[inode].displacement[1],
                         boundary_nodes[inode].displacement[2]);
  }
  /* Calculate the weights with the solver that fits the support of the RBF. */
  const Eigen::MatrixXd alpha = rbf_function->is_compactly_supported ()
                                  ? t8_rbf_solve_compactly_supported (boundary_nodes, *rbf_function, displacements)
                                  : t8_rbf_solve_globally_supported (boundary_nodes, *rbf_function, displacements);
  /* Copy the calculated weights back to the boundary nodes. */
  for (size_t inode = 0; inode < num_boundary_nodes; ++inode) {
    boundary_nodes[inode].weight[0] = alpha (inode, 0);
    boundary_nodes[inode].weight[1] = alpha (inode, 1);
    boundary_nodes[inode].weight[2] = alpha (inode, 2);
  }
}
