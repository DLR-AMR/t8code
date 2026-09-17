#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <span>
#include <vector>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/io/vtk_shape.hxx"

namespace t8_mra
{

namespace detail
{

/// Barycentric index of the node at a given position in VTK's Lagrange triangle
/// ordering (port of vtkHigherOrderTriangle::BarycentricIndex).
[[nodiscard]] inline std::array<int, 3>
vtk_triangle_barycentric_index (int index, int order)
{
  int max = order;
  int min = 0;

  /// Scope into the correct inner triangle
  while (index != 0 && index >= 3 * order) {
    index -= 3 * order;
    max -= 2;
    ++min;
    order -= 3;
  }

  if (index == 0)
    return { min, min, max };
  if (index == 1)
    return { max, min, min };
  if (index == 2)
    return { min, max, min };

  std::array<int, 3> bindex;
  index -= 3;
  const int dim = index / (order - 1);
  const int offset = index - dim * (order - 1);
  bindex[dim] = min + 1 + offset;
  bindex[(dim + 1) % 3] = min;
  bindex[(dim + 2) % 3] = max - 1 - offset;

  return bindex;
}

}  // namespace detail

/// Triangle: barycentric Lagrange node layout and the affine map over the three
/// corners. The corner order is data-driven, since Bey refinement permutes it.
template <>
struct vtk_shape<T8_ECLASS_TRIANGLE>
{
  static constexpr int DIM = 2;
  static constexpr int MAX_LAGRANGE_ORDER = 10;
  using point = std::array<double, DIM>;

  [[nodiscard]] static std::vector<point>
  lagrange_nodes (int order)
  {
    std::vector<point> nodes (shape_traits<T8_ECLASS_TRIANGLE>::dof (order + 1));

    for (auto i = 0u; i < nodes.size (); ++i) {
      const auto bindex = detail::vtk_triangle_barycentric_index (static_cast<int> (i), order);
      nodes[i] = { static_cast<double> (bindex[0]) / order, static_cast<double> (bindex[1]) / order };
    }

    return nodes;
  }

  [[nodiscard]] static int
  vertex_slot (int corner, const std::array<int, 3> &order)
  {
    return order[corner];
  }

  /** @brief Affine map in barycentric weights over the VTK-ordered corners. */
  [[nodiscard]] static std::array<double, 3>
  to_physical (const point &ref, std::span<const std::array<double, 3>> vertices)
  {
    const double w0 = 1.0 - ref[0] - ref[1];
    const double w1 = ref[0];
    const double w2 = ref[1];

    std::array<double, 3> phys {};
    for (auto d = 0; d < 3; ++d)
      phys[d] = w0 * vertices[0][d] + w1 * vertices[1][d] + w2 * vertices[2][d];

    return phys;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
