#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>
#include <cmath>

///TODO Maybe nested namespace for cartesian
namespace t8_mra
{

/**
 * @brief Maps a point from physical interval [xL, xR] to reference interval [0, 1]
 *
 * Performs affine transformation: x_ref = (x - xL) / (xR - xL)
 *
 * @param x Point in physical interval [xL, xR]
 * @param xL Left endpoint of physical interval
 * @param xR Right endpoint of physical interval
 * @return double Point in reference interval [0, 1]
 */
inline double
ref_1d (double x, double xL, double xR)
{
  return (x - xL) / (xR - xL);
}

/**
 * @brief Maps a point from reference interval [0, 1] to physical interval [xL, xR]
 *
 * Performs affine transformation: x = x_ref * (xR - xL) + xL
 *
 * @param x_ref Point in reference interval [0, 1]
 * @param xL Left endpoint of physical interval
 * @param xR Right endpoint of physical interval
 * @return double Point in physical interval [xL, xR]
 */
inline double
deref_1d (double x_ref, double xL, double xR)
{
  return x_ref * (xR - xL) + xL;
}

/**
 * @brief Maps a multi-dimensional point from physical element to reference element [0,1]^DIM
 *
 * For cartesian elements, this is a dimension-wise affine transformation.
 *
 * @tparam DIM Spatial dimension
 * @param x Point in physical element
 * @param vertices_min Lower corner of physical element (coordinates of vertex 0)
 * @param vertices_max Upper corner of physical element (coordinates of opposite vertex)
 * @return std::array<double, DIM> Point in reference element [0,1]^DIM
 */
template <unsigned int DIM>
inline std::array<double, DIM>
ref (const std::array<double, DIM> &x, const std::array<double, DIM> &vertices_min,
     const std::array<double, DIM> &vertices_max)
{
  std::array<double, DIM> x_ref;
  for (unsigned int d = 0; d < DIM; ++d)
    x_ref[d] = ref_1d (x[d], vertices_min[d], vertices_max[d]);

  return x_ref;
}

/**
 * @brief Maps a multi-dimensional point from reference element [0,1]^DIM to physical element
 *
 * For cartesian elements, this is a dimension-wise affine transformation.
 *
 * @tparam DIM Spatial dimension
 * @param x_ref Point in reference element [0,1]^DIM
 * @param vertices_min Lower corner of physical element
 * @param vertices_max Upper corner of physical element
 * @return std::array<double, DIM> Point in physical element
 */
template <unsigned int DIM>
inline std::array<double, DIM>
deref (const std::array<double, DIM> &x_ref, const std::array<double, DIM> &vertices_min,
       const std::array<double, DIM> &vertices_max)
{
  std::array<double, DIM> x;
  for (unsigned int d = 0; d < DIM; ++d)
    x[d] = deref_1d (x_ref[d], vertices_min[d], vertices_max[d]);

  return x;
}

/**
 * @brief Computes the Jacobian determinant of the mapping from reference to physical element
 *
 * For cartesian elements with axis-aligned edges, the Jacobian is diagonal:
 * J = diag(h_0, h_1, ..., h_{DIM-1}) where h_i = vertices_max[i] - vertices_min[i]
 * det(J) = product of edge lengths
 *
 * @tparam DIM Spatial dimension
 * @param vertices_min Lower corner of physical element
 * @param vertices_max Upper corner of physical element
 * @return double Absolute value of Jacobian determinant
 */
template <unsigned int DIM>
inline double
jacobian_determinant (const std::array<double, DIM> &vertices_min, const std::array<double, DIM> &vertices_max)
{
  double det = 1.0;
  for (unsigned int d = 0; d < DIM; ++d)
    det *= (vertices_max[d] - vertices_min[d]);

  return std::abs (det);
}

/**
 * @brief Extracts vertices from t8code vertex array for cartesian elements
 *
 * For cartesian elements (LINE, QUAD, HEX), extracts the min and max vertices
 * which define the axis-aligned bounding box.
 *
 * @tparam DIM Spatial dimension
 * @param physical_vertices Array of vertex coordinates from t8code (layout depends on element type)
 * @param vertices_min Output: lower corner coordinates
 * @param vertices_max Output: upper corner coordinates
 */
template <unsigned int DIM>
inline void
extract_cartesian_vertices (const double physical_vertices[][3], std::array<double, DIM> &vertices_min,
                            std::array<double, DIM> &vertices_max)
{
  if constexpr (DIM == 1) {
    // LINE: vertices[0] and vertices[1]
    vertices_min[0] = physical_vertices[0][0];
    vertices_max[0] = physical_vertices[1][0];
  }
  else if constexpr (DIM == 2) {
    // QUAD: t8code uses row-major ordering:
    // vertex 0: (xmin, ymin), vertex 1: (xmax, ymin)
    // vertex 2: (xmin, ymax), vertex 3: (xmax, ymax)
    vertices_min[0] = physical_vertices[0][0];
    vertices_min[1] = physical_vertices[0][1];
    vertices_max[0] = physical_vertices[3][0];
    vertices_max[1] = physical_vertices[3][1];
  }
  else if constexpr (DIM == 3) {
    // HEX: vertex 0 is (xmin, ymin, zmin), vertex 6 is (xmax, ymax, zmax)
    vertices_min[0] = physical_vertices[0][0];
    vertices_min[1] = physical_vertices[0][1];
    vertices_min[2] = physical_vertices[0][2];
    vertices_max[0] = physical_vertices[6][0];
    vertices_max[1] = physical_vertices[6][1];
    vertices_max[2] = physical_vertices[6][2];
  }
}

/**
 * @brief Transforms quadrature points from reference to physical element
 *
 * Takes quadrature points on [0,1]^DIM and maps them to the physical element.
 *
 * @tparam DIM Spatial dimension
 * @param ref_quad_points Reference quadrature points (flattened: [x0,y0,z0, x1,y1,z1, ...])
 * @param num_points Number of quadrature points
 * @param vertices_min Lower corner of physical element
 * @param vertices_max Upper corner of physical element
 * @return std::vector<double> Physical quadrature points (flattened, same layout)
 */
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

    std::array<double, DIM> x_phys = deref<DIM> (x_ref, vertices_min, vertices_max);

    for (unsigned int d = 0; d < DIM; ++d)
      phys_quad_points[DIM * i + d] = x_phys[d];
  }

  return phys_quad_points;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
