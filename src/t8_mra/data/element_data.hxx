#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include "t8_eclass/t8_eclass.h"

#include "t8_mra/data/levelindex_map.hxx"
#include "t8_mra/data/levelmultiindex.hxx"

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

template <t8_eclass TShape>
static constexpr unsigned short
t8_eclass_dim ()
{
  if constexpr (TShape == T8_ECLASS_LINE)
    return 1;
  else if constexpr (TShape == T8_ECLASS_TRIANGLE || TShape == T8_ECLASS_QUAD)
    return 2;
  else if constexpr (TShape == T8_ECLASS_TET || TShape == T8_ECLASS_HEX || TShape == T8_ECLASS_PRISM
                     || TShape == T8_ECLASS_PYRAMID)
    return 3;
  else
    static_assert (false, "Invalid element class");
}

template <t8_eclass TShape>
static constexpr unsigned short
t8_eclass_num_children ()
{
  if constexpr (TShape == T8_ECLASS_LINE)
    return 2;
  else if constexpr (TShape == T8_ECLASS_TRIANGLE || TShape == T8_ECLASS_QUAD)
    return 4;
  else if constexpr (TShape == T8_ECLASS_TET || TShape == T8_ECLASS_HEX || TShape == T8_ECLASS_PRISM
                     || TShape == T8_ECLASS_PYRAMID)
    return 8;  /// TODO
  else
    static_assert (false, "Invalid element class");
}

/**
 * @brief Compute number of DOF for a given element type and polynomial order
 *
 * For simplex elements (TRIANGLE, TET): uses binomial coefficient
 * For cartesian elements (LINE, QUAD, HEX): uses tensor product P^DIM
 */
template <t8_eclass TShape, unsigned short P>
static constexpr unsigned short
compute_dof ()
{
  constexpr unsigned short DIM = t8_eclass_dim<TShape> ();

  if constexpr (TShape == T8_ECLASS_TRIANGLE || TShape == T8_ECLASS_TET)
    return binom (DIM + P - 1, DIM);
  else if constexpr (TShape == T8_ECLASS_LINE)
    return P;
  else if constexpr (TShape == T8_ECLASS_QUAD)
    return P * P;
  else if constexpr (TShape == T8_ECLASS_HEX)
    return P * P * P;
  else {
    static_assert (false, "Unsupported element type for DOF computation");
  }
}

/// TODO template specialization
template <t8_eclass TShape, unsigned short U, unsigned short P>
struct element_data
{

  static constexpr t8_eclass Shape = TShape;
  static constexpr unsigned short DIM = t8_eclass_dim<TShape> ();
  static constexpr unsigned short NUM_CHILDREN = t8_eclass_num_children<TShape> ();
  static constexpr unsigned short U_DIM = U;

  static constexpr unsigned short P_DIM = P;
  static constexpr unsigned short DOF = compute_dof<TShape, P> ();
  static constexpr unsigned short W_DOF = DOF * NUM_CHILDREN;

  // Fixed-size storage: keeps the struct trivially copyable, so element
  // data can be shipped between ranks as raw bytes (repartitioning, ghost
  // exchange) and lives inline in the dense maps.
  std::array<double, U_DIM * DOF> u_coeffs = {};  // Single-scale coefficients
  /// TODO get rid of details
  std::array<double, U_DIM * W_DOF> d_coeffs = {};  // Detail coefficients
  double vol = 0.0;

  std::array<int, 3> order = {};  // Point order

  size_t static dg_idx (size_t u, size_t p) noexcept
  {
    return u * DOF + p;
  }

  size_t static wavelet_idx (size_t k, size_t u, size_t p) noexcept
  {
    return k * U_DIM * DOF + u * DOF + p;
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
