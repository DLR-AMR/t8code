#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include "t8_eclass/t8_eclass.h"

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

/// Per-cell DG data. Shape facts come from shape_traits<TShape>.
template <t8_eclass TShape, unsigned short U, unsigned short P>
struct element_data
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr unsigned short DIM = shape_traits<TShape>::DIM;
  static constexpr unsigned short NUM_CHILDREN = shape_traits<TShape>::NUM_CHILDREN;
  static constexpr unsigned short U_DIM = U;

  static constexpr unsigned short P_DIM = P;
  static constexpr unsigned short DOF = shape_traits<TShape>::dof (P);
  static constexpr unsigned short W_DOF = DOF * NUM_CHILDREN;

  // Fixed-size storage: keeps the struct trivially copyable, so element
  // data can be shipped between ranks as raw bytes (repartitioning, ghost
  // exchange).
  std::array<double, U_DIM * DOF> u_coeffs = {};  // Single-scale coefficients
  double vol = 0.0;

  std::array<int, 3> order = {};  // Point order

  size_t static dg_idx (size_t u, size_t p) noexcept
  {
    return u * DOF + p;
  }
};

/// Leaf data plus detail (wavelet) coefficients.
template <t8_eclass TShape, unsigned short U, unsigned short P>
struct detail_data: element_data<TShape, U, P>
{
  using base = element_data<TShape, U, P>;

  std::array<double, base::U_DIM * base::W_DOF> d_coeffs = {};

  size_t static wavelet_idx (size_t k, size_t u, size_t p) noexcept
  {
    return k * base::U_DIM * base::DOF + u * base::DOF + p;
  }
};

template <typename T>
struct forest_data
{
  using lmi_type = levelmultiindex<T::Shape>;

  sc_array_t *lmi_idx;
  t8_mra::levelindex_map<lmi_type, T> *lmi_map;

  void *mra_instance;  // Pointer to multiscale object for callbacks
};

template <typename T>
t8_mra::levelmultiindex<T::Shape>
get_lmi_from_forest_data (const t8_mra::forest_data<T> *forest_data, size_t idx)
{
  return *reinterpret_cast<t8_mra::levelmultiindex<T::Shape> *> (t8_sc_array_index_locidx (forest_data->lmi_idx, idx));
}

template <typename T>
void
set_lmi_forest_data (t8_mra::forest_data<T> *forest_data, size_t idx, const t8_mra::levelmultiindex<T::Shape> &lmi)
{
  *reinterpret_cast<t8_mra::levelmultiindex<T::Shape> *> (t8_sc_array_index_locidx (forest_data->lmi_idx, idx)) = lmi;
}

}  // namespace t8_mra

#endif
