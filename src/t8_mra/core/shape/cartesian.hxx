#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_eclass/t8_eclass.h"

#include "t8_mra/core/shape/mst_policy.hxx"
#include "t8_mra/core/shape_traits.hxx"

namespace t8_mra
{

template <>
struct shape_traits<T8_ECLASS_LINE>
{
  static constexpr unsigned short DIM = 1;
  static constexpr unsigned short NUM_CHILDREN = 2;
  static constexpr int NUM_VERTICES = 2;
  static constexpr int VTK_CELL_TYPE = 68;  // VTK_LAGRANGE_CURVE

  [[nodiscard]] static constexpr unsigned short
  dof (unsigned short P)
  {
    return P;
  }
};

template <>
struct shape_traits<T8_ECLASS_QUAD>
{
  static constexpr unsigned short DIM = 2;
  static constexpr unsigned short NUM_CHILDREN = 4;
  static constexpr int NUM_VERTICES = 4;
  static constexpr int VTK_CELL_TYPE = 70;  // VTK_LAGRANGE_QUADRILATERAL

  [[nodiscard]] static constexpr unsigned short
  dof (unsigned short P)
  {
    return P * P;
  }
};

template <>
struct shape_traits<T8_ECLASS_HEX>
{
  static constexpr unsigned short DIM = 3;
  static constexpr unsigned short NUM_CHILDREN = 8;
  static constexpr int NUM_VERTICES = 8;
  static constexpr int VTK_CELL_TYPE = 72;  // VTK_LAGRANGE_HEXAHEDRON

  [[nodiscard]] static constexpr unsigned short
  dof (unsigned short P)
  {
    return P * P * P;
  }
};

/// Cartesian shapes carry no vertex-order information.
template <t8_eclass TShape>
  requires is_cartesian<TShape>
struct ordering_policy<TShape>
{
  template <typename TData>
  static void
  adjust_parent_order (TData & /*unused*/)
  {
  }

  template <typename TData>
  static void
  adjust_child_order (TData & /*unused*/, int /*unused*/, const TData & /*unused*/)
  {
  }
};

/// Cartesian MST scaling: L2-orthonormal reference basis
template <t8_eclass TShape>
  requires is_cartesian<TShape>
struct mst_scaling_policy<TShape>
{
  [[nodiscard]] static constexpr double
  forward_scaling_factor (unsigned int num_children)
  {
    return 1.0 / static_cast<double> (num_children);
  }

  [[nodiscard]] static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }

  [[nodiscard]] static constexpr double
  detail_norm_scale (double /*unused*/)
  {
    return 1.0;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
