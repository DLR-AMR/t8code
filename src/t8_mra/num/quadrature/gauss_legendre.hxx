#pragma once

#ifdef T8_ENABLE_MRA

#include <gsl/gsl_integration.h>
#include <memory>
#include <vector>
#include <stdexcept>

namespace t8_mra
{

/**
 * @brief Generates 1D Gauss-Legendre quadrature points and weights on [0,1].
 *
 * The rule integrates polynomials up to degree 2*num_points-1 exactly.
 * points and weights are overwritten with num_points entries.
 */
inline void
gauss_legendre_1d (int num_points, std::vector<double> &points, std::vector<double> &weights)
{
  if (num_points <= 0)
    throw std::invalid_argument ("Number of quadrature points must be positive");

  const std::unique_ptr<gsl_integration_fixed_workspace, decltype (&gsl_integration_fixed_free)> workspace (
    gsl_integration_fixed_alloc (gsl_integration_fixed_legendre, num_points, 0.0, 1.0, 0.0, 0.0),
    gsl_integration_fixed_free);

  if (!workspace)
    throw std::runtime_error ("Failed to allocate GSL integration workspace");

  const double *nodes = gsl_integration_fixed_nodes (workspace.get ());
  const double *w = gsl_integration_fixed_weights (workspace.get ());

  points.assign (nodes, nodes + num_points);
  weights.assign (w, w + num_points);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
