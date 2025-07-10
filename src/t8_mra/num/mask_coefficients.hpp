#pragma once

#ifdef T8_ENABLE_MRA

#include <vector>
#include <stdexcept>

#include <t8_mra/num/mat.hpp>
#include <t8_eclass.h>

namespace t8_mra
{
template <t8_eclass TShape>
void
initialize (size_t polynomial_degree, std::vector<t8_mra::mat>& mask_coeffs, std::vector<t8_mra::mat>& inv_mask_coeffs)
{
  throw std::out_of_range ("Element shape is not supported in "
                           "t8_mra::mask_coefficients::initialize");
}

template <>
void
initialize<T8_ECLASS_TRIANGLE> (size_t polynomial_degree, std::vector<t8_mra::mat>& mask_coeffs,
                                std::vector<t8_mra::mat>& inv_mask_coeffs);

}  // namespace t8_mra

#include <t8_mra/num/mask_coeffs_triangle.hpp>

#endif
