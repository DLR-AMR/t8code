#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <span>
#include <vector>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/basis.hxx"
#include "t8_mra/num/mat.hxx"

namespace t8_mra
{

/**
 * @brief Converts nodal DG data to the modal coefficients MRA stores.
 *
 * Given a nodal value n_j = u(x_j) at reference nodes x_j, the modal
 * coefficients solve V m = n with the Vandermonde V_ji = phi_i(x_j). The
 * cartesian basis is orthonormal with unit normalization on every cell, so V is
 * cell-independent: it is built and LU-factored once at construction, then each
 * per-cell conversion is a pair of triangular solves per component.
 *
 * Nodes must be DOF distinct reference points making V nonsingular (e.g. a
 * tensor Gauss-Lobatto / equispaced nodal set of order P). Nodal and modal
 * buffers are component-major: index u*DOF + j, matching element_data u_coeffs.
 */
template <t8_eclass TShape, unsigned int U, unsigned int P>
  requires is_cartesian<TShape>
class nodal_to_modal {
 public:
  static constexpr unsigned int DIM = shape_traits<TShape>::DIM;
  static constexpr unsigned int DOF = shape_traits<TShape>::dof (P);

  explicit nodal_to_modal (const std::array<std::array<double, DIM>, DOF> &nodes): vandermonde (DOF, DOF), perm (DOF)
  {
    for (auto j = 0u; j < DOF; ++j) {
      const auto phi = basis<TShape, P>::eval (nodes[j]);
      for (auto i = 0u; i < DOF; ++i)
        vandermonde (j, i) = phi[i];
    }
    lu_factors (vandermonde, perm);
  }

  /// Convert one cell's nodal values (U*DOF) into modal coeffs (U*DOF).
  void
  operator() (std::span<const double> nodal, std::span<double> modal) const
  {
    std::array<double, DOF> rhs;
    for (auto u = 0u; u < U; ++u) {
      const auto off = u * DOF;
      for (auto j = 0u; j < DOF; ++j)
        rhs[j] = nodal[off + j];
      lu_solve (vandermonde, perm, rhs);
      for (auto i = 0u; i < DOF; ++i)
        modal[off + i] = rhs[i];
    }
  }

  std::array<double, U * DOF>
  operator() (std::span<const double> nodal) const
  {
    std::array<double, U * DOF> modal;
    (*this) (nodal, modal);

    return modal;
  }

 private:
  mat vandermonde;           // LU-factored Vandermonde phi_i(x_j)
  std::vector<size_t> perm;  // pivot permutation
};

/**
 * @brief Reconstructs nodal DG values from modal coefficients (inverse of
 * nodal_to_modal).
 *
 * n_j = sum_i m_i phi_i(x_j), a matvec against the Vandermonde built from the
 * reference `nodes`; phi_i(x_j) is evaluated once at construction. Cartesian
 * only, buffers component-major (index u*DOF + j).
 */
template <t8_eclass TShape, unsigned int U, unsigned int P>
  requires is_cartesian<TShape>
class modal_to_nodal {
 public:
  static constexpr unsigned int DIM = shape_traits<TShape>::DIM;
  static constexpr unsigned int DOF = shape_traits<TShape>::dof (P);

  explicit modal_to_nodal (const std::array<std::array<double, DIM>, DOF> &nodes)
  {
    for (auto j = 0u; j < DOF; ++j)
      phi_at_node[j] = basis<TShape, P>::eval (nodes[j]);
  }

  void
  operator() (std::span<const double> modal, std::span<double> nodal) const
  {
    for (auto u = 0u; u < U; ++u) {
      const auto off = u * DOF;
      for (auto j = 0u; j < DOF; ++j) {
        double v = 0.0;
        for (auto i = 0u; i < DOF; ++i)
          v += modal[off + i] * phi_at_node[j][i];
        nodal[off + j] = v;
      }
    }
  }

  std::array<double, U * DOF>
  operator() (std::span<const double> modal) const
  {
    std::array<double, U * DOF> nodal;
    (*this) (modal, nodal);

    return nodal;
  }

 private:
  std::array<std::array<double, DOF>, DOF> phi_at_node;  // phi_i(x_j)
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
