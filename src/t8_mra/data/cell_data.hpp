#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <array>
#include "t8_eclass.h"

#include "t8_mra/data/levelindex_map.hpp"
#include "t8_mra/data/levelmultiindex.hpp"

namespace t8_mra
{

constexpr inline size_t
binom (int n, int k) noexcept
{
  return (k > n)                  ? 0
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
  static constexpr t8_eclass Shape = TShape;
  static constexpr unsigned short DIM = 2;  /// TODO
  static constexpr unsigned short U_DIM = U;

  static constexpr unsigned short P_DIM = P;
  static constexpr unsigned short DOF = binom (DIM + P_DIM - 1, DIM);
  static constexpr unsigned short W_DOF = DOF * 3;

  std::vector<double> u_coeffs;  // Single-scale coefficients
  std::vector<double> d_coeffs;  // Detail coefficients
  bool significant;
  std::array<int, 3> order;  // Point order

  explicit data_per_element ()
    : u_coeffs (U_DIM * DOF, {}), d_coeffs (U_DIM * W_DOF, {}), significant (false), order ({})
  {
  }
};

// template <t8_eclass TShape>
// struct element_data
// {
//   t8_mra::levelmultiindex<TShape> lmi_idx;
// };

template <typename T>
struct forest_data
{
  sc_array_t *lmi_idx;
  t8_mra::levelindex_map<T> *lmi_map;
};

template <typename T>
t8_mra::levelmultiindex<T::Shape>
get_lmi_from_forest_data (const t8_mra::forest_data<T> *forest_data, size_t idx)
{
  return *reinterpret_cast<t8_mra::levelmultiindex<T::Shape> *> (t8_sc_array_index_locidx (forest_data->lmi_idx, idx));
}

}  // namespace t8_mra

#endif
