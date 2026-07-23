/**
 * @file
 * @brief Dunavant quadrature rules over the reference triangle.
 * Coefficients from David Dunavant, "High Degree Efficient Symmetrical
 * Gaussian Quadrature Rules for the Triangle", IJNME 21 (1985), 1129-1148.
 * Original C tables by John Burkardt; the data lives in dunavant_table.hxx.
 */

#pragma once

#ifdef T8_ENABLE_MRA

#include <span>
#include <vector>

namespace t8_mra
{

/// A reference-triangle quadrature rule: points flattened as [x0,y0, x1,y1, ...]
/// with one weight per point.
struct dunavant_quadrature
{
  std::vector<double> points;
  std::vector<double> weights;
};

/// Expand Dunavant rule `rule` (1..20) into its full point/weight set.
dunavant_quadrature
dunavant_rule (int rule);

/// Map reference-triangle points to physical space. tri holds the three
/// vertices [x0,y0, x1,y1, x2,y2]; ref holds points [x0,y0, ...]; the returned
/// vector holds the physical points in the same flattened layout.
std::vector<double>
reference_to_physical_t3 (std::span<const double> tri, std::span<const double> ref);

}  // namespace t8_mra

#endif
