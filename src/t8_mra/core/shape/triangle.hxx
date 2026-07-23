#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/core/shape/mst_policy.hxx"
#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/data/triangle_order.hxx"

namespace t8_mra
{

template <>
struct shape_traits<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned short DIM = 2;
  static constexpr unsigned short NUM_CHILDREN = 4;
  static constexpr int NUM_VERTICES = 3;
  static constexpr int VTK_CELL_TYPE = 69;  // VTK_LAGRANGE_TRIANGLE

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return binom (DIM + P - 1, DIM);
  }
};

/// Triangle vertex order follows the Bey refinement type across levels.
template <>
struct ordering_policy<T8_ECLASS_TRIANGLE>
{
  template <typename T>
  static void
  adjust_parent_order (T &data)
  {
    triangle_order::get_parent_order (data.order);
  }

  template <typename T>
  static void
  adjust_child_order (T &child_data, int child_id, const T &parent_data)
  {
    child_data.order = parent_data.order;
    triangle_order::get_point_order (child_data.order, child_id);
  }
};

/// Triangle MST scaling: no child averaging in the forward transform; the
/// detail norm scales by 1/vol.
template <>
struct mst_scaling_policy<T8_ECLASS_TRIANGLE>
{
  static constexpr double
  forward_scaling_factor (unsigned int /*num_children*/)
  {
    return 1.0;
  }

  static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }

  static constexpr double
  detail_norm_scale (double vol)
  {
    return 1.0 / vol;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
