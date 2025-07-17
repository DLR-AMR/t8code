#pragma once

#ifdef T8_ENABLE_MRA

#include <array>

#include <t8.h>
#include <t8_forest/t8_forest_general.h>

namespace t8_mra
{

template <t8_eclass TShape, int U>
struct quadrature
{
  int order_num;
  template <typename TFunc>
  std::array<double, U>
  projection (t8_forest_t forest, int tree_idx, const t8_element_t* elem, TFunc&& func)
  {
  }
};

}  // namespace t8_mra

#endif
