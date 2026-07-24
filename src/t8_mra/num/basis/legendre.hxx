#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <cmath>

#include <gsl/gsl_sf_legendre.h>

namespace t8_mra
{

/**
 * @brief The p-th Legendre polynomial at x, shifted to [0,1] and L2-normalized.
 *
 * Shift x_std = 2x-1 and multiply by sqrt(2p+1) for orthonormality on [0,1].
 *
 * @param x Point in [0,1].
 * @param p Polynomial degree.
 */
[[nodiscard]] inline double
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
 * @brief Derivative of the p-th normalized Legendre polynomial at x on [0,1].
 *
 * @tparam P number of 1D modes (degrees 0..P-1); sizes the GSL scratch exactly.
 * @param p Polynomial degree in [0, P).
 *
 * GSL writes degrees 0..p; the chain-rule factor 2 (from x_std = 2x-1) and the
 * sqrt(2p+1) normalization are applied to the requested degree.
 */
template <int P>
[[nodiscard]] inline double
phi_prime_1d (double x, int p)
{
  const auto x_std = 2.0 * x - 1.0;

  // Derivative of the constant mode is zero (and keeps p+1 <= P for p > 0).
  if (p == 0)
    return 0.0;

  std::array<double, P> leg_array;
  std::array<double, P> deriv_array;
  gsl_sf_legendre_Pl_deriv_array (p, x_std, leg_array.data (), deriv_array.data ());

  return 2.0 * deriv_array[p] * std::sqrt (2.0 * p + 1.0);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
