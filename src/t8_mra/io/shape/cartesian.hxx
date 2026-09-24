#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <span>
#include <vector>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/io/vtk_shape.hxx"

namespace t8_mra
{

namespace detail
{

/// Index of the grid node (i, j) in VTK's Lagrange quad ordering
/// (port of vtkHigherOrderQuadrilateral::PointIndexFromIJK, uniform order).
[[nodiscard]] inline int
vtk_quad_point_index (int i, int j, int order)
{
  const bool ibdy = (i == 0 || i == order);
  const bool jbdy = (j == 0 || j == order);

  if (ibdy && jbdy)  /// vertex
    return (i != 0 ? (j != 0 ? 2 : 1) : (j != 0 ? 3 : 0));

  int offset = 4;
  if (ibdy || jbdy) {  /// edge
    if (!ibdy)
      return (i - 1) + ((j != 0) ? 2 * (order - 1) : 0) + offset;
    return (j - 1) + ((i != 0) ? order - 1 : 3 * (order - 1)) + offset;
  }

  offset += 4 * (order - 1);  /// interior
  return offset + (i - 1) + (order - 1) * (j - 1);
}

/// Index of the grid node (i, j, k) in VTK's Lagrange hex ordering
/// (port of vtkHigherOrderHexahedron::PointIndexFromIJK, uniform order).
[[nodiscard]] inline int
vtk_hex_point_index (int i, int j, int k, int order)
{
  const bool ibdy = (i == 0 || i == order);
  const bool jbdy = (j == 0 || j == order);
  const bool kbdy = (k == 0 || k == order);
  const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3)  /// vertex
    return ((i != 0) ? ((j != 0) ? 2 : 1) : ((j != 0) ? 3 : 0)) + ((k != 0) ? 4 : 0);

  int offset = 8;
  if (nbdy == 2) {  /// edge
    if (!ibdy)
      return (i - 1) + ((j != 0) ? 2 * (order - 1) : 0) + ((k != 0) ? 4 * (order - 1) : 0) + offset;
    if (!jbdy)
      return (j - 1) + ((i != 0) ? order - 1 : 3 * (order - 1)) + ((k != 0) ? 4 * (order - 1) : 0) + offset;
    offset += 8 * (order - 1);
    return (k - 1) + (order - 1) * ((i != 0) ? ((j != 0) ? 2 : 1) : ((j != 0) ? 3 : 0)) + offset;
  }

  offset += 12 * (order - 1);
  const int face_size = (order - 1) * (order - 1);
  if (nbdy == 1) {  /// face
    if (ibdy)
      return (j - 1) + (order - 1) * (k - 1) + ((i != 0) ? face_size : 0) + offset;
    offset += 2 * face_size;
    if (jbdy)
      return (i - 1) + (order - 1) * (k - 1) + ((j != 0) ? face_size : 0) + offset;
    offset += 2 * face_size;
    return (i - 1) + (order - 1) * (j - 1) + ((k != 0) ? face_size : 0) + offset;
  }

  offset += 6 * face_size;  /// interior
  return offset + (i - 1) + (order - 1) * ((j - 1) + (order - 1) * (k - 1));
}

/// Inverse of a corner permutation: slot of each t8code corner.
template <size_t N>
[[nodiscard]] constexpr std::array<int, N>
inverse_permutation (const std::array<int, N> &perm) noexcept
{
  std::array<int, N> inverse {};
  for (auto slot = 0u; slot < N; ++slot)
    inverse[perm[slot]] = static_cast<int> (slot);

  return inverse;
}

/// t8code corner -> VTK slot for a cartesian shape.
/// t8code numbers corners with the first coordinate fastest; VTK walks the face cycle.
template <t8_eclass TShape>
[[nodiscard]] constexpr auto
vtk_corner_permutation ()
{
  if constexpr (TShape == T8_ECLASS_LINE)
    return std::array<int, 2> { 0, 1 };
  else if constexpr (TShape == T8_ECLASS_QUAD)
    return std::array<int, 4> { 0, 1, 3, 2 };
  else  /// HEX
    return std::array<int, 8> { 0, 1, 3, 2, 4, 5, 7, 6 };
}

/// Reference coords of each VTK-ordered corner, for the multilinear map.
template <t8_eclass TShape>
[[nodiscard]] constexpr auto
vtk_corner_coords ()
{
  if constexpr (TShape == T8_ECLASS_LINE)
    return std::array<std::array<double, 1>, 2> { { { 0 }, { 1 } } };
  else if constexpr (TShape == T8_ECLASS_QUAD)
    return std::array<std::array<double, 2>, 4> { { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } } };
  else  /// HEX
    return std::array<std::array<double, 3>, 8> {
      { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } }
    };
}

}  // namespace detail

/// Cartesian shapes (LINE, QUAD, HEX): tensor Lagrange node layout and the
/// multilinear geometry map over the corners.
template <t8_eclass TShape>
  requires is_cartesian<TShape>
struct vtk_shape<TShape>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;
  static constexpr int MAX_LAGRANGE_ORDER = 10;
  using point = std::array<double, DIM>;

  [[nodiscard]] static std::vector<point>
  lagrange_nodes (int order)
  {
    std::vector<point> nodes (shape_traits<TShape>::dof (order + 1));

    if constexpr (DIM == 1) {
      nodes[0] = { 0.0 };
      nodes[1] = { 1.0 };
      for (int i = 1; i < order; ++i)
        nodes[i + 1] = { static_cast<double> (i) / order };
    }
    else if constexpr (DIM == 2) {
      for (int j = 0; j <= order; ++j)
        for (int i = 0; i <= order; ++i)
          nodes[detail::vtk_quad_point_index (i, j, order)]
            = { static_cast<double> (i) / order, static_cast<double> (j) / order };
    }
    else {
      for (int k = 0; k <= order; ++k)
        for (int j = 0; j <= order; ++j)
          for (int i = 0; i <= order; ++i)
            nodes[detail::vtk_hex_point_index (i, j, k, order)]
              = { static_cast<double> (i) / order, static_cast<double> (j) / order, static_cast<double> (k) / order };
    }

    return nodes;
  }

  /// Inverse of the corner permutation, so the caller can scatter by t8code corner.
  [[nodiscard]] static int
  vertex_slot (int corner, const std::array<int, 3> & /*unused*/)
  {
    static constexpr auto slot = detail::inverse_permutation (detail::vtk_corner_permutation<TShape> ());

    return slot[corner];
  }

  /** @brief Multilinear map: tensor of (1-x_d)/x_d weights per VTK-ordered corner. */
  [[nodiscard]] static std::array<double, 3>
  to_physical (const point &ref, std::span<const std::array<double, 3>> vertices)
  {
    constexpr auto corners = detail::vtk_corner_coords<TShape> ();
    std::array<double, 3> phys = { 0.0, 0.0, 0.0 };

    for (auto v = 0u; v < corners.size (); ++v) {
      double weight = 1.0;
      for (auto d = 0; d < DIM; ++d)
        weight *= corners[v][d] == 1.0 ? ref[d] : 1.0 - ref[d];

      for (auto d = 0; d < 3; ++d)
        phys[d] += weight * vertices[v][d];
    }

    return phys;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
