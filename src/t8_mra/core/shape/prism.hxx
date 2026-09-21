#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_eclass/t8_eclass.h"

#include "t8_mra/core/shape/mst_policy.hxx"
#include "t8_mra/core/shape/triangle.hxx"
#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/data/triangle_order.hxx"

namespace t8_mra
{

template <>
struct shape_traits<T8_ECLASS_PRISM>
{
  static constexpr unsigned short DIM = 3;
  static constexpr unsigned short NUM_CHILDREN = 8;
  static constexpr int NUM_VERTICES = 6;
  static constexpr int VTK_CELL_TYPE = 73;  /// VTK_LAGRANGE_WEDGE

  [[nodiscard]] static constexpr unsigned short
  dof (unsigned short P)
  {
    return shape_traits<T8_ECLASS_TRIANGLE>::dof (P) * P;
  }
};

/// Prism children factor as triangle x line (t8_dprism_child), so only the
/// triangle half carries vertex order and the line half contributes nothing.
template <>
struct ordering_policy<T8_ECLASS_PRISM>
{
  template <typename TData>
  static void
  adjust_parent_order (TData &data)
  {
    triangle_order::get_parent_order (data.order);
  }

  template <typename TData>
  static void
  adjust_child_order (TData &child_data, int child_id, const TData &parent_data)
  {
    child_data.order = parent_data.order;

    triangle_order::get_point_order (child_data.order, child_id % shape_traits<T8_ECLASS_TRIANGLE>::NUM_CHILDREN);
  }
};

template <>
struct mst_scaling_policy<T8_ECLASS_PRISM>
{
  [[nodiscard]] static constexpr double
  forward_scaling_factor (unsigned int /*unused*/)
  {
    return 1.0;
  }

  [[nodiscard]] static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }

  [[nodiscard]] static constexpr double
  detail_norm_scale (double vol)
  {
    return 1.0 / vol;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
