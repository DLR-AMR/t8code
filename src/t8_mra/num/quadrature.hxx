#pragma once

#ifdef T8_ENABLE_MRA

#include <gsl/gsl_integration.h>
#include <vector>
#include <stdexcept>

namespace t8_mra
{

/**
 * @brief Generates 1D Gauss-Legendre quadrature points and weights on [0,1]
 *
 * Uses GSL to generate Gauss-Legendre quadrature nodes and weights.
 * The quadrature rule integrates polynomials up to degree 2*n-1 exactly,
 * where n is the number of quadrature points.
 *
 * @param num_points Number of quadrature points
 * @param points Output vector for quadrature points (will be resized)
 * @param weights Output vector for quadrature weights (will be resized)
 */
inline void
gauss_legendre_1d (int num_points, std::vector<double> &points, std::vector<double> &weights)
{
  if (num_points <= 0)
    throw std::invalid_argument ("Number of quadrature points must be positive");

  // Allocate GSL workspace for Gauss-Legendre quadrature on [0,1]
  gsl_integration_fixed_workspace *workspace
    = gsl_integration_fixed_alloc (gsl_integration_fixed_legendre, num_points, 0.0, 1.0, 0.0, 0.0);

  if (!workspace)
    throw std::runtime_error ("Failed to allocate GSL integration workspace");

  // Get pointers to nodes and weights
  double *nodes_ptr = gsl_integration_fixed_nodes (workspace);
  double *weights_ptr = gsl_integration_fixed_weights (workspace);

  // Copy to output vectors
  points.resize (num_points);
  weights.resize (num_points);

  for (int i = 0; i < num_points; ++i) {
    points[i] = nodes_ptr[i];
    weights[i] = weights_ptr[i];
  }

  // Free workspace
  gsl_integration_fixed_free (workspace);
}

/**
 * @brief Generates 1D Gauss-Legendre quadrature on arbitrary interval [a,b]
 *
 * @param num_points Number of quadrature points
 * @param a Left endpoint of interval
 * @param b Right endpoint of interval
 * @param points Output vector for quadrature points (will be resized)
 * @param weights Output vector for quadrature weights (will be resized)
 */
inline void
gauss_legendre_1d_interval (int num_points, double a, double b, std::vector<double> &points,
                            std::vector<double> &weights)
{
  if (a >= b)
    throw std::invalid_argument ("Interval endpoints must satisfy a < b");

  // Allocate GSL workspace for Gauss-Legendre quadrature on [a,b]
  gsl_integration_fixed_workspace *workspace
    = gsl_integration_fixed_alloc (gsl_integration_fixed_legendre, num_points, a, b, 0.0, 0.0);

  if (!workspace)
    throw std::runtime_error ("Failed to allocate GSL integration workspace");

  // Get pointers to nodes and weights
  double *nodes_ptr = gsl_integration_fixed_nodes (workspace);
  double *weights_ptr = gsl_integration_fixed_weights (workspace);

  // Copy to output vectors
  points.resize (num_points);
  weights.resize (num_points);

  for (int i = 0; i < num_points; ++i) {
    points[i] = nodes_ptr[i];
    weights[i] = weights_ptr[i];
  }

  // Free workspace
  gsl_integration_fixed_free (workspace);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
