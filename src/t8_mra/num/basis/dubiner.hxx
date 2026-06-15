#pragma once
#ifdef T8_ENABLE_MRA

#include <array>
#include <cmath>

namespace t8_mra
{
/// Jacobi polynomial P_n^{(alpha,beta)}(x) by three-term recurrence.
double
jacobi (int n, double alpha, double beta, double x);

/// Derivative d/dx P_n^{(alpha,beta)}(x) = (n+alpha+beta+1)/2 P_{n-1}^{(alpha+1,beta+1)}(x).
inline double
jacobi_deriv (int n, double alpha, double beta, double x)
{
  if (n == 0)
    return 0.0;
  return 0.5 * (n + alpha + beta + 1.0) * jacobi (n - 1, alpha + 1.0, beta + 1.0, x);
}

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

/// Reference gradient {d/dtau1, d/dtau2} of the I-th Dubiner scaling function.
/// Evaluated at interior (Dunavant) points; the collapsed-coordinate gradient
/// is singular only at the apex (c -> 0, p == 1), guarded by clamping c.
template <int I>
std::array<double, 2>
scaling_function_gradient (double tau1, double tau2)
{
  constexpr int d = detail::dubiner_degree (I);
  constexpr int p = I - d * (d + 1) / 2;
  constexpr int q = d - p;

  const double norm = std::sqrt (2.0 * (2 * p + 1) * (p + q + 1));
  const double c = 1.0 - tau2;
  const double cc = c < 1e-12 ? 1e-12 : c;
  const double a = 1.0 - 2.0 * tau1 / cc;
  const double b = 2.0 * tau2 - 1.0;

  const double pa = jacobi (p, 0.0, 0.0, a);
  const double pb = jacobi (q, 2.0 * p + 1.0, 0.0, b);
  const double pad = jacobi_deriv (p, 0.0, 0.0, a);
  const double pbd = jacobi_deriv (q, 2.0 * p + 1.0, 0.0, b);

  double dt1 = 0.0;
  double dt2 = norm * 2.0 * pa * pbd * detail::power<p> (cc);  // term from d(P_q(b))/dtau2

  // The d/dtau1 path and the c-power derivatives only contribute for p >= 1
  // (for p == 0 the Jacobi-in-a derivative pad is zero anyway).
  if constexpr (p >= 1) {
    const double cpm1 = detail::power<p - 1> (cc);
    dt1 = norm * pad * (-2.0) * cpm1 * pb;        // d a/d tau1 = -2/c
    dt2 += norm * (-static_cast<double> (p)) * pa * pb * cpm1;  // d(c^p)/d tau2
    if constexpr (p >= 2)
      dt2 += norm * (-2.0) * tau1 * pad * pb * detail::power<p - 2> (cc);  // d a/d tau2 = -2 tau1/c^2
    else
      dt2 += norm * (-2.0) * tau1 * pad * pb / cc;  // p == 1: c^{p-2} = 1/c
  }

  return { dt1, dt2 };
}

}  // namespace t8_mra
#endif
