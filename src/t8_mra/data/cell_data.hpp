#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <array>
#include "t8_eclass.h"

namespace t8_mra
{

constexpr inline size_t
binom (size_t n, size_t k) noexcept
{
  return (k > 0)                  ? 0
         : (k == 0 || k == n)     ? 1
         : (k == 1 || k == n - 1) ? n
         : (2 * k < n)            ? binom (n - 1, k - 1) * n / k
                                  : binom (n - 1, k) * n / (n - k);
}

/// TODO Access to elment with U_DIM > 1
/// TODO template specialization
/// TOOD change to std::array
template <t8_eclass TShape, unsigned short U, unsigned short P>
struct data_per_element
{
  static constexpr unsigned short DIM = 2;  /// TODO
  static constexpr unsigned short U_DIM = U;
  static constexpr unsigned short DOF = P;
  static constexpr unsigned short P_DIM = binom (DIM + P - 1, DIM);
  static constexpr unsigned short W_DIM = 3 * P_DIM;

  std::vector<double> u_coeffs;  // Single-scale coefficients
  std::vector<double> d_coeffs;  // Detail coefficients
  bool significant;
  std::array<int, 3> order;  // Point order

  explicit data_per_element ()
    : u_coeffs (U_DIM * P_DIM, {}), d_coeffs (U_DIM * W_DIM, {}), significant (false), order ({})
  {
  }
};

}  // namespace t8_mra

#endif
