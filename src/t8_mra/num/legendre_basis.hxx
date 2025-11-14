#pragma once

#ifdef T8_ENABLE_MRA

#include <gsl/gsl_sf_legendre.h>
#include <vector>
#include <cmath>

namespace t8_mra
{

/**
 * @brief Evaluates the p-th Legendre polynomial at point x on [0,1]
 *
 * The Legendre polynomials are shifted from the standard [-1,1] interval
 * to [0,1] and L2-normalized.
 *
 * For orthonormality on [0,1]:
 * - Legendre on [-1,1] with norm sqrt(2/(2p+1))
 * - Shift to [0,1]: x_std = 2x-1, dx = d(x_std)/2
 * - This gives normalization factor sqrt(2(2p+1)) for [0,1]
 *
 * @param x Point in [0,1] where to evaluate
 * @param p Polynomial degree (0, 1, 2, ...)
 * @return double Value of the p-th normalized Legendre polynomial at x
 */
inline double
phi_1d (double x, int p)
{
  // Transform x from [0,1] to [-1,1] for standard Legendre polynomials
  double x_std = 2.0 * x - 1.0;

  // Evaluate Legendre polynomial using GSL
  double leg_value = gsl_sf_legendre_Pl (p, x_std);

  // Apply L2-normalization for [0,1] interval: sqrt(2*p + 1)
  return leg_value * std::sqrt (2.0 * p + 1.0);
}

/**
 * @brief Evaluates the derivative of the p-th Legendre polynomial at point x on [0,1]
 *
 * Computes d/dx of the p-th Legendre polynomial, taking into account the
 * transformation from [0,1] to [-1,1] and the L2-normalization.
 *
 * @param x Point in [0,1] where to evaluate the derivative
 * @param p Polynomial degree (0, 1, 2, ...)
 * @return double Value of the derivative at x
 */
inline double
phi_prime_1d (double x, int p)
{
  // Transform x from [0,1] to [-1,1]
  const auto x_std = 2.0 * x - 1.0;

  // Derivative of constant polynomial is zero
  if (p == 0)
    return 0.0;

  std::vector<double> leg_array (p + 1);
  std::vector<double> deriv_array (p + 1);

  // Compute all Legendre polynomials and derivatives up to order p
  gsl_sf_legendre_Pl_deriv_array (p, x_std, leg_array.data (), deriv_array.data ());
  const auto deriv_std = deriv_array[p];

  // Apply chain rule: d/dx_01 = d/dx_std * dx_std/dx_01
  // dx_std/dx_01 = 2.0 (from x_std = 2*x - 1)
  // Also apply L2-normalization factor sqrt(2*p + 1) for [0,1]
  return 2.0 * deriv_std * std::sqrt (2.0 * p + 1.0);
}

/**
 * @brief Evaluates all Legendre basis functions up to order p_max at point x
 *
 * @param x Point in [0,1] where to evaluate
 * @param p_max Maximum polynomial degree
 * @param values Output array of size (p_max + 1) to store the values
 */
inline void
phi_1d_array (double x, int p_max, double *values)
{
  // Transform x from [0,1] to [-1,1]
  const auto x_std = 2.0 * x - 1.0;

  // Compute all Legendre polynomials up to p_max
  std::vector<double> leg_array (p_max + 1);
  gsl_sf_legendre_Pl_array (p_max, x_std, leg_array.data ());

  // Apply L2-normalization for [0,1] interval to each polynomial
  for (int p = 0; p <= p_max; ++p)
    values[p] = leg_array[p] * std::sqrt (2.0 * p + 1.0);
}

/**
 * @brief Evaluates derivatives of all Legendre basis functions up to order p_max at point x
 *
 * @param x Point in [0,1] where to evaluate
 * @param p_max Maximum polynomial degree
 * @param derivs Output array of size (p_max + 1) to store the derivatives
 */
inline void
phi_prime_1d_array (double x, int p_max, double *derivs)
{
  // Transform x from [0,1] to [-1,1]
  const auto x_std = 2.0 * x - 1.0;

  // Compute all Legendre polynomials and derivatives up to p_max
  std::vector<double> leg_array (p_max + 1);
  std::vector<double> deriv_array (p_max + 1);
  gsl_sf_legendre_Pl_deriv_array (p_max, x_std, leg_array.data (), deriv_array.data ());

  // Apply chain rule and L2-normalization for [0,1] interval to each derivative
  for (int p = 0; p <= p_max; ++p)
    derivs[p] = 2.0 * deriv_array[p] * std::sqrt (2.0 * p + 1.0);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
