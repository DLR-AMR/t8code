#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <cmath>

namespace t8_mra
{

namespace detail
{

/// Unnormalized Legendre polynomials of degrees 0..P-1 at xi in [-1,1].
template <int P>
[[nodiscard]] inline std::array<double, P>
legendre_raw (double xi) noexcept
{
  std::array<double, P> value {};
  value[0] = 1.0;

  if constexpr (P > 1)
    value[1] = xi;

  for (auto p = 1; p + 1 < P; ++p)
    value[p + 1] = ((2 * p + 1) * xi * value[p] - p * value[p - 1]) / (p + 1);

  return value;
}

/// sqrt(2p+1), the L2 normalization of degree p on [0,1].
template <int P>
inline const std::array<double, P> legendre_norm = [] {
  std::array<double, P> norm {};
  for (auto p = 0; p < P; ++p)
    norm[p] = std::sqrt (2.0 * p + 1.0);

  return norm;
}();

}  // namespace detail

/// Normalized Legendre values and derivatives of degrees 0..P-1 at one point.
template <int P>
struct legendre_modes
{
  std::array<double, P> value {};
  std::array<double, P> derivative {};
};

/** @brief Degrees 0..P-1 at x in [0,1], orthonormal there, in one recurrence. */
template <int P>
[[nodiscard]] inline std::array<double, P>
legendre_values (double x) noexcept
{
  auto value = detail::legendre_raw<P> (2.0 * x - 1.0);

  for (auto p = 0; p < P; ++p)
    value[p] *= detail::legendre_norm<P>[p];

  return value;
}

/** @brief Degrees 0..P-1 and their derivatives at x in [0,1], orthonormal there. */
template <int P>
[[nodiscard]] inline legendre_modes<P>
legendre_at (double x) noexcept
{
  const auto raw = detail::legendre_raw<P> (2.0 * x - 1.0);
  legendre_modes<P> modes;

  /// P'_{p+1} = P'_{p-1} + (2p+1) P_p on [-1,1]; xi = 2x-1
  if constexpr (P > 1)
    modes.derivative[1] = 1.0;

  for (auto p = 1; p + 1 < P; ++p)
    modes.derivative[p + 1] = modes.derivative[p - 1] + (2 * p + 1) * raw[p];

  for (auto p = 0; p < P; ++p) {
    modes.value[p] = raw[p] * detail::legendre_norm<P>[p];
    modes.derivative[p] *= 2.0 * detail::legendre_norm<P>[p];
  }

  return modes;
}

/** @brief The p-th Legendre polynomial at x, shifted to [0,1] and L2-normalized. */
[[nodiscard]] inline double
phi_1d (double x, int p) noexcept
{
  if (p == 0)
    return 1.0;

  const auto xi = 2.0 * x - 1.0;
  auto previous = 1.0;
  auto current = xi;

  for (auto n = 1; n < p; ++n) {
    const auto next = ((2 * n + 1) * xi * current - n * previous) / (n + 1);
    previous = current;
    current = next;
  }

  return current * std::sqrt (2.0 * p + 1.0);
}

/** @brief Derivative of the p-th normalized Legendre polynomial at x on [0,1]. */
[[nodiscard]] inline double
phi_prime_1d (double x, int p) noexcept
{
  if (p == 0)
    return 0.0;

  const auto xi = 2.0 * x - 1.0;
  auto previous = 1.0;
  auto current = xi;
  auto derivative_previous = 0.0;
  auto derivative = 1.0;

  for (auto n = 1; n < p; ++n) {
    const auto next = ((2 * n + 1) * xi * current - n * previous) / (n + 1);
    const auto next_derivative = derivative_previous + (2 * n + 1) * current;

    previous = current;
    current = next;
    derivative_previous = derivative;
    derivative = next_derivative;
  }

  return 2.0 * derivative * std::sqrt (2.0 * p + 1.0);
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
