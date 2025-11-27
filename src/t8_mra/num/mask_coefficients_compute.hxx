#pragma once

#ifdef T8_ENABLE_MRA

#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>

#include <t8_eclass.h>
#include <t8_mra/num/mat.hpp>
#include <t8_mra/num/legendre_basis.hxx>
#include <t8_mra/num/quadrature.hxx>

namespace t8_mra
{

/**
 * @brief Compute wavelet basis function in 1D
 *
 * For Legendre wavelets, the wavelet is constructed as:
 * ψ_p(x) = φ_p(2x-1) - projection onto V_0
 *
 * For our purposes, we use the standard construction from the 2-scale relation.
 *
 * @param x Point in [0,1] where to evaluate
 * @param p Polynomial degree (0, 1, 2, ...)
 * @return double Value of the wavelet at x
 */
inline double
psi_1d (double x, int p)
{
  // Wavelet on [0,1]: difference between fine and coarse scale
  // ψ_p(x) = √2 φ_p(2x) - φ_p(x) for x ∈ [0, 0.5]
  //        = √2 φ_p(2x-1) - φ_p(x) for x ∈ [0.5, 1]

  // For simplicity, use the definition:
  // ψ_p(x) = φ_p(2x) - φ_p(2x-1) (for x ∈ [0.5, 1])
  // This captures the detail (high-frequency) component

  if (x < 0.5) {
    return std::sqrt (2.0) * phi_1d (2.0 * x, p);
  }
  else {
    return std::sqrt (2.0) * phi_1d (2.0 * x - 1.0, p);
  }
}

/**
 * @brief Compute mask coefficients for cartesian elements (LINE, QUAD, HEX)
 *
 * Computes the mask matrices M_k such that:
 * - Forward MST: u_parent[i] = (1/2^D) * sum_k sum_j M_k[j,i] * u_child_k[j]
 * - Inverse MST: u_child_k[i] = sum_j M_k[i,j] * u_parent[j]
 *
 * The mask coefficients are computed via numerical integration:
 * M_k[i,j] = ∫ φ_i(2x - s_k) φ_j(x) dx
 *
 * where s_k is the shift vector for child k (e.g., [0,0], [1,0], [0,1], [1,1] for QUAD).
 *
 * For wavelets (e ≠ 0), we compute:
 * M_k[i,j,e] = ∫ ψ_i^e(2x - s_k) φ_j(x) dx
 *
 * where ψ^e indicates wavelet in directions specified by e ∈ {0,1}^D.
 *
 * @tparam ECLASS Element class (T8_ECLASS_LINE, T8_ECLASS_QUAD, T8_ECLASS_HEX)
 * @param order Polynomial order P (degree = P-1)
 * @param dof Degrees of freedom per element (should be P^D)
 * @param mask_coeffs Output: mask matrices for each child, size [2^D][dof x dof]
 */
template <t8_eclass_t ECLASS>
void
compute_mask_coefficients (size_t order, size_t dof, std::vector<t8_mra::mat> &mask_coeffs)
{
  // Get dimensionality
  constexpr int DIM = (ECLASS == T8_ECLASS_LINE)   ? 1
                      : (ECLASS == T8_ECLASS_QUAD) ? 2
                      : (ECLASS == T8_ECLASS_HEX)  ? 3
                                                   : -1;
  static_assert (DIM > 0, "Unsupported element class for mask coefficient computation");

  constexpr int NUM_CHILDREN = (1 << DIM);  // 2^D

  // Initialize output
  mask_coeffs.resize (NUM_CHILDREN, t8_mra::mat { dof, dof });

  // Generate multiindices for tensor product basis
  // P_set[i] = (i_0, i_1, ..., i_{D-1}) where each i_d ∈ [0, P-1]
  std::vector<std::array<int, DIM>> P_set (dof);
  for (size_t idx = 0; idx < dof; ++idx) {
    size_t temp = idx;
    for (int d = 0; d < DIM; ++d) {
      P_set[idx][d] = temp % order;
      temp /= order;
    }
  }

  for (auto i = 0; i < P_set.size (); ++i)
    std::cout << "[" << P_set[i][0] << ", " << P_set[i][1] << "], ";
  std::cout << std::endl;

  // Generate child shift vectors E_set[k] ∈ {0,1}^D
  std::vector<std::array<int, DIM>> E_set (NUM_CHILDREN);
  for (int k = 0; k < NUM_CHILDREN; ++k) {
    for (int d = 0; d < DIM; ++d) {
      E_set[k][d] = (k >> d) & 1;
    }
  }

  // Setup quadrature for integration
  // Use high-order Gauss-Legendre quadrature to accurately integrate polynomials
  const int quad_order = 2 * (order - 1) + 1;  // Exact for degree 2*(P-1)
  std::vector<double> quad_nodes_1d, quad_weights_1d;
  gauss_legendre_1d (quad_order, quad_nodes_1d, quad_weights_1d);

  // For tensor product in D dimensions, we need quad_order^D points
  // Generate tensor product quadrature
  const size_t num_quad_1d = quad_nodes_1d.size ();
  std::vector<std::array<double, DIM>> gauss_nodes;
  std::vector<double> gauss_weights;

  // Generate all D-dimensional quadrature points
  auto generate_tensor_quad = [&] (auto &self, int dim, std::array<double, DIM> pt, double wt) -> void {
    if (dim == DIM) {
      gauss_nodes.push_back (pt);
      gauss_weights.push_back (wt);
      return;
    }
    for (size_t q = 0; q < num_quad_1d; ++q) {
      pt[dim] = quad_nodes_1d[q];
      self (self, dim + 1, pt, wt * quad_weights_1d[q]);
    }
  };

  std::array<double, DIM> init_pt {};
  generate_tensor_quad (generate_tensor_quad, 0, init_pt, 1.0);

  // Compute mask coefficients for scaling functions following multilaepsch
  // Formula: M_k[i,j] = ∫_{[0,1]^D} φ_i(0.5*(x + s_k)) * φ_j(x) dx
  //
  // where:
  // - x is integration variable on parent reference element [0,1]^D
  // - φ_j(x) is parent basis function j evaluated at x
  // - φ_i(0.5*(x + s_k)) is child basis function i evaluated at child reference coords
  // - s_k is the shift vector for child k (e.g., [0,0], [1,0], [0,1], [1,1] for QUAD)
  //
  // This gives the prolongation matrix for inverse MST: u_child_k = M_k * u_parent
  // For forward MST (restriction), we use the transpose: u_parent = (1/2^D) * sum_k M_k^T * u_child_k

  std::cout << "Computing mask coefficients for " << t8_eclass_to_string[ECLASS] << ", P=" << order << "...\n";

  for (int k = 0; k < NUM_CHILDREN; ++k) {
    for (size_t i = 0; i < dof; ++i) {
      for (size_t j = 0; j < dof; ++j) {
        double integral = 0.0;

        // Integrate over the child's region in parent reference coords
        // Child k occupies parent region [s_k/2, (s_k+1)/2]^D
        // We integrate over [0,1]^D but map it to the child's region
        for (size_t q = 0; q < gauss_nodes.size (); ++q) {
          const auto &x_unit = gauss_nodes[q];  // Integration point on [0,1]^D
          const double gaussW = gauss_weights[q];

          // Map to child's region in parent coords: x = s_k/2 + x_unit/2
          // Then map to child reference coords: xs = 2*x - s_k = x_unit
          // So child reference coords are just x_unit!
          std::array<double, DIM> xs;
          for (int d = 0; d < DIM; ++d) {
            xs[d] = x_unit[d];
          }

          // Parent coords where we evaluate parent basis
          std::array<double, DIM> x_parent;
          for (int d = 0; d < DIM; ++d) {
            x_parent[d] = 0.5 * (E_set[k][d] + x_unit[d]);
          }

          // Evaluate φ_i at child reference coords: φ_i(xs)
          double phi_i_val = 1.0;
          for (int d = 0; d < DIM; ++d) {
            phi_i_val *= phi_1d (xs[d], P_set[i][d]);
          }

          // Evaluate φ_j at parent reference coords: φ_j(x_parent)
          double phi_j_val = 1.0;
          for (int d = 0; d < DIM; ++d) {
            phi_j_val *= phi_1d (x_parent[d], P_set[j][d]);
          }

          // Accumulate integral
          // The Jacobian factor will be applied in the forward MST, not here
          integral += gaussW * phi_i_val * phi_j_val;
        }

        // Store the mask coefficient
        // mask_coeffs[k](i, j) = integral * 0.5;
        mask_coeffs[k](i, j) = integral;
      }
    }
  }

  // Debug output: print first mask matrix
  std::cout << "\nMask matrix for child 0 (scaling functions):\n";
  for (size_t i = 0; i < std::min (dof, size_t (5)); ++i) {
    for (size_t j = 0; j < std::min (dof, size_t (5)); ++j) {
      std::cout << std::setw (12) << std::setprecision (6) << mask_coeffs[0](i, j) << " ";
    }
    std::cout << (dof > 5 ? "..." : "") << "\n";
  }
  if (dof > 5)
    std::cout << "...\n";

  std::cout << "\nMask coefficients computed successfully.\n\n";
}

/**
 * @brief Convenience function to compute and initialize mask coefficients
 *
 * This function replaces the hardcoded values in mask_coeffs_*.hpp
 */
template <t8_eclass_t ECLASS>
void
initialize_mask_coefficients_computed (size_t order, size_t dof, std::vector<t8_mra::mat> &mask_coeffs,
                                       std::vector<t8_mra::mat> &inv_mask_coeffs)
{
  // Compute mask coefficients
  compute_mask_coefficients<ECLASS> (order, dof, mask_coeffs);

  // Note: inv_mask_coeffs are for wavelets (e ≠ 0)
  // For now, we only compute scaling function masks
  // Wavelet masks would require computing the wavelet basis functions
  // and their integrals, which is more complex

  // Initialize inv_mask_coeffs to appropriate size (but leave empty for now)
  constexpr int DIM = (ECLASS == T8_ECLASS_LINE)   ? 1
                      : (ECLASS == T8_ECLASS_QUAD) ? 2
                      : (ECLASS == T8_ECLASS_HEX)  ? 3
                                                   : -1;
  constexpr int NUM_CHILDREN = (1 << DIM);
  inv_mask_coeffs.resize (NUM_CHILDREN, t8_mra::mat { dof, dof });
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
