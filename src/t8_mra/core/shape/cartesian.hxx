#pragma once

#ifdef T8_ENABLE_MRA

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

  static constexpr unsigned short
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

  static constexpr unsigned short
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

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return P * P * P;
  }
};

/// Cartesian shapes carry no vertex-order information: no-op.
template <t8_eclass Shape>
  requires is_cartesian<Shape>
struct ordering_policy<Shape>
{
  template <typename T>
  static void
  adjust_parent_order (T &)
  {
  }

  template <typename T>
  static void
  adjust_child_order (T &, int, const T &)
  {
  }
};

/// Cartesian MST scaling: L2-orthonormal reference basis, so the forward
/// transform averages over children and the detail norm needs no scaling.
template <t8_eclass Shape>
  requires is_cartesian<Shape>
struct mst_scaling_policy<Shape>
{
  static constexpr double
  forward_scaling_factor (unsigned int num_children)
  {
    return 1.0 / static_cast<double> (num_children);
  }

  static constexpr double
  inverse_scaling_factor ()
  {
    return 1.0;
  }

  static constexpr double
  detail_norm_scale (double /*vol*/)
  {
    return 1.0;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
