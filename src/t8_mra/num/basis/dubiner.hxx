#pragma once
#ifdef T8_ENABLE_MRA

#include <cmath>

namespace t8_mra
{
/// Jacobi polynomial P_n^{(alpha,beta)}(x) by three-term recurrence.
double
jacobi (int n, double alpha, double beta, double x);

namespace detail
{
/// Total degree d of linear basis index i (smallest d with (d+1)(d+2)/2 > i).
constexpr int
dubiner_degree (int i)
{
  int d = 0;
  while ((d + 1) * (d + 2) / 2 <= i)
    ++d;
  return d;
}

/// c^N for small compile-time N
template <int N>
constexpr double
power (double c)
{
  double r = 1.0;
  for (int k = 0; k < N; ++k)
    r *= c;
  return r;
}
}  // namespace detail

/// I-th orthonormal Dubiner scaling function on the reference triangle (area
/// 1/2), via Jacobi recurrence; valid for any I. tau1, tau2 in [0,1].
template <int I>
double
scaling_function (double tau1, double tau2)
{
  constexpr int d = detail::dubiner_degree (I);
  constexpr int p = I - d * (d + 1) / 2;
  constexpr int q = d - p;

  // Collapsed coordinates. At the apex (tau2 -> 1) the c^p factor kills every
  // p>=1 term.
  const double c = 1.0 - tau2;
  const double a = c < 1e-12 ? -1.0 : 1.0 - 2.0 * tau1 / c;
  const double b = 2.0 * tau2 - 1.0;
  const double norm = std::sqrt (2.0 * (2 * p + 1) * (p + q + 1));

  return norm * jacobi (p, 0.0, 0.0, a) * detail::power<p> (c) * jacobi (q, 2.0 * p + 1.0, 0.0, b);
}

}  // namespace t8_mra
#endif
