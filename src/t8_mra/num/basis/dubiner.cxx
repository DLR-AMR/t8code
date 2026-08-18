#ifdef T8_ENABLE_MRA

#include "t8_mra/num/basis/dubiner.hxx"

namespace t8_mra
{
double
jacobi (int n, double alpha, double beta, double x)
{
  if (n == 0)
    return 1.0;
  double p_prev = 1.0;
  double p = 0.5 * (alpha - beta) + 0.5 * (alpha + beta + 2.0) * x;
  for (int k = 2; k <= n; ++k) {
    const double c = 2.0 * k + alpha + beta;
    const double a1 = 2.0 * k * (k + alpha + beta) * (c - 2.0);
    const double a2 = (c - 1.0) * (alpha * alpha - beta * beta);
    const double a3 = (c - 1.0) * c * (c - 2.0);
    const double a4 = 2.0 * (k + alpha - 1.0) * (k + beta - 1.0) * c;
    const double p_next = ((a2 + a3 * x) * p - a4 * p_prev) / a1;
    p_prev = p;
    p = p_next;
  }
  return p;
}

}  // namespace t8_mra

#endif
