#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>

namespace t8_mra
{

/// Affine map from reference [0,1] to physical [xL, xR].
constexpr double
deref_1d (double x_ref, double xL, double xR) noexcept
{
  return x_ref * (xR - xL) + xL;
}

/// Dimension-wise affine map [0,1]^DIM -> physical cartesian cell.
template <unsigned int DIM>
constexpr std::array<double, DIM>
deref (const std::array<double, DIM> &x_ref, const std::array<double, DIM> &vertices_min,
       const std::array<double, DIM> &vertices_max) noexcept
{
  std::array<double, DIM> x;
  for (unsigned int d = 0; d < DIM; ++d)
    x[d] = deref_1d (x_ref[d], vertices_min[d], vertices_max[d]);
  return x;
}

/// Min/max corners of an axis-aligned cartesian cell from t8code vertices. The
/// vertex permutation applied in dg/cartesian.hxx puts the lower corner at
/// index 0 and the upper corner at the last index (2/3/7 for LINE/QUAD/HEX).
template <unsigned int DIM>
inline void
extract_cartesian_vertices (const double physical_vertices[][3], std::array<double, DIM> &vertices_min,
                            std::array<double, DIM> &vertices_max) noexcept
{
  constexpr int max_vertex = (DIM == 1) ? 1 : (DIM == 2) ? 2 : 7;
  for (unsigned int d = 0; d < DIM; ++d) {
    vertices_min[d] = physical_vertices[0][d];
    vertices_max[d] = physical_vertices[max_vertex][d];
  }
}

/// Maps the flattened reference quadrature points ([x0,y0,..., x1,y1,...]) to
/// the physical cartesian cell.
template <unsigned int DIM>
inline std::vector<double>
transform_quad_points (const std::vector<double> &ref_quad_points, size_t num_points,
                       const std::array<double, DIM> &vertices_min, const std::array<double, DIM> &vertices_max)
{
  std::vector<double> phys_quad_points (DIM * num_points);

  for (size_t i = 0; i < num_points; ++i) {
    std::array<double, DIM> x_ref;
    for (unsigned int d = 0; d < DIM; ++d)
      x_ref[d] = ref_quad_points[DIM * i + d];

    const std::array<double, DIM> x_phys = deref<DIM> (x_ref, vertices_min, vertices_max);

    for (unsigned int d = 0; d < DIM; ++d)
      phys_quad_points[DIM * i + d] = x_phys[d];
  }

  return phys_quad_points;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
