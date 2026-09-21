#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <span>
#include <vector>

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/io/shape/triangle.hxx"
#include "t8_mra/io/vtk_shape.hxx"

namespace t8_mra
{

namespace detail
{

/// Reference coords of the six wedge corners: triangle corner (vertex % 3) at
/// height (vertex / 3), matching t8_dprism_vertex_integer_coords.
inline constexpr std::array<std::array<double, 3>, 6> wedge_corner
  = { { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 } } };

/// vtkWedge edge list: bottom triangle, top triangle, then the three verticals.
inline constexpr std::array<std::array<int, 2>, 9> wedge_edge
  = { { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 3, 4 }, { 4, 5 }, { 5, 3 }, { 0, 3 }, { 1, 4 }, { 2, 5 } } };

/// vtkWedge triangle faces. The top face is listed 3,5,4, so its local axes are
/// swapped against the bottom one.
inline constexpr std::array<std::array<int, 3>, 2> wedge_tri_face = { { { 0, 1, 2 }, { 3, 5, 4 } } };

/// vtkWedge quad faces. Local i runs c0->c1 (vertical), local j runs c0->c3.
inline constexpr std::array<std::array<int, 4>, 3> wedge_quad_face
  = { { { 0, 3, 4, 1 }, { 1, 4, 5, 2 }, { 2, 5, 3, 0 } } };

/// Reference point of a triangle-face node, from its barycentric index. The
/// weights are cyclic (b2 on c0, b0 on c1, b1 on c2), as in
/// vtk_triangle_barycentric_index.
[[nodiscard]] inline std::array<double, 3>
wedge_tri_face_point (const std::array<int, 3> &face, const std::array<int, 3> &bindex, int order)
{
  std::array<double, 3> p {};
  for (auto d = 0; d < 3; ++d)
    p[d] = (bindex[2] * wedge_corner[face[0]][d] + bindex[0] * wedge_corner[face[1]][d]
            + bindex[1] * wedge_corner[face[2]][d])
           / order;

  return p;
}

}  // namespace detail

/// Prism: VTK's Lagrange wedge node layout (6 corners, 9 edges, 5 faces,
/// interior) and the geometry map, barycentric in the base and linear in z.
template <>
struct vtk_shape<T8_ECLASS_PRISM>
{
  static constexpr int DIM = 3;
  static constexpr int MAX_LAGRANGE_ORDER = 10;
  using point = std::array<double, DIM>;

  [[nodiscard]] static std::vector<point>
  lagrange_nodes (int order)
  {
    const int tri_interior = (order - 1) * (order - 2) / 2;
    std::vector<point> nodes;
    nodes.reserve (shape_traits<T8_ECLASS_PRISM>::dof (order + 1));

    for (const auto &corner : detail::wedge_corner)
      nodes.push_back (corner);

    for (const auto &edge : detail::wedge_edge) {
      const auto &from = detail::wedge_corner[edge[0]];
      const auto &to = detail::wedge_corner[edge[1]];

      for (auto i = 1; i < order; ++i) {
        const double s = static_cast<double> (i) / order;
        nodes.push_back (
          { from[0] + s * (to[0] - from[0]), from[1] + s * (to[1] - from[1]), from[2] + s * (to[2] - from[2]) });
      }
    }

    for (const auto &face : detail::wedge_tri_face)
      for (auto m = 0; m < tri_interior; ++m)
        nodes.push_back (
          detail::wedge_tri_face_point (face, detail::vtk_triangle_barycentric_index (3 * order + m, order), order));

    // Quad-face interior in the face's own ordering: local i fastest.
    for (const auto &face : detail::wedge_quad_face)
      for (auto j = 1; j < order; ++j)
        for (auto i = 1; i < order; ++i) {
          const auto u = static_cast<double> (i) / order;
          const auto v = static_cast<double> (j) / order;
          std::array<double, 3> p {};

          for (auto d = 0; d < 3; ++d)
            p[d] = (1.0 - u) * (1.0 - v) * detail::wedge_corner[face[0]][d]
                   + u * (1.0 - v) * detail::wedge_corner[face[1]][d] + u * v * detail::wedge_corner[face[2]][d]
                   + (1.0 - u) * v * detail::wedge_corner[face[3]][d];

          nodes.push_back (p);
        }

    // Cell interior: triangle-interior index fastest, then the layer.
    for (auto k = 1; k < order; ++k)
      for (auto m = 0; m < tri_interior; ++m) {
        const auto bindex = detail::vtk_triangle_barycentric_index (3 * order + m, order);
        nodes.push_back ({ static_cast<double> (bindex[0]) / order, static_cast<double> (bindex[1]) / order,
                           static_cast<double> (k) / order });
      }

    T8_ASSERT (nodes.size () == shape_traits<T8_ECLASS_PRISM>::dof (order + 1));

    return nodes;
  }

  /// Both triangle ends share the order, and t8code pairs them by +3.
  [[nodiscard]] static int
  vertex_slot (int corner, const std::array<int, 3> &order)
  {
    return order[corner % 3] + 3 * (corner / 3);
  }

  /** @brief Barycentric in the base triangle, linear between the two ends. */
  [[nodiscard]] static std::array<double, 3>
  to_physical (const point &ref, std::span<const std::array<double, 3>> vertices)
  {
    const auto w0 = 1.0 - ref[0] - ref[1];
    const auto w1 = ref[0];
    const auto w2 = ref[1];
    const auto z = ref[2];

    std::array<double, 3> phys {};
    for (auto d = 0; d < 3; ++d)
      phys[d] = (1.0 - z) * (w0 * vertices[0][d] + w1 * vertices[1][d] + w2 * vertices[2][d])
                + z * (w0 * vertices[3][d] + w1 * vertices[4][d] + w2 * vertices[5][d]);

    return phys;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
