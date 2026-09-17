/**
 * @file
 * @brief Dunavant quadrature rules over the reference triangle.
 * Coefficients from David Dunavant, "High Degree Efficient Symmetrical
 * Gaussian Quadrature Rules for the Triangle", IJNME 21 (1985), 1129-1148.
 * Original C tables by John Burkardt; the data lives in dunavant_table.hxx.
 */

#pragma once

#ifdef T8_ENABLE_MRA

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "t8_mra/num/quadrature/dunavant_table.hxx"

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
[[nodiscard]] inline dunavant_quadrature
dunavant_rule (int rule)
{
  if (rule < 1 || rule > 20)
    throw std::out_of_range ("dunavant_rule: rule must be in [1, 20]");

  dunavant_quadrature quad;

  const auto add = [&quad] (double x, double y, double weight) {
    quad.points.push_back (x);
    quad.points.push_back (y);
    quad.weights.push_back (weight);
  };

  // Each orbit expands into mult points by cyclic permutation of its
  // barycentric coordinates (x = bary[k], y = bary[k+1], third coord implied).
  for (const auto &orbit : dunavant_table::rule (rule)) {
    const auto &bary = orbit.bary;
    switch (orbit.mult) {
    case 1:
      add (bary[0], bary[1], orbit.weight);
      break;
    case 3:
      for (int k = 0; k < 3; ++k)
        add (bary[k], bary[(k + 1) % 3], orbit.weight);
      break;
    case 6:
      for (int k = 0; k < 3; ++k)
        add (bary[k], bary[(k + 1) % 3], orbit.weight);
      for (int k = 0; k < 3; ++k)
        add (bary[(k + 1) % 3], bary[k], orbit.weight);
      break;
    default:
      throw std::logic_error ("dunavant_rule: invalid orbit multiplicity");
    }
  }

  return quad;
}

/// Map reference-triangle points to physical space. tri holds the three
/// vertices [x0,y0, x1,y1, x2,y2]; ref holds points [x0,y0, ...]; the returned
/// vector holds the physical points in the same flattened layout.
[[nodiscard]] inline std::vector<double>
reference_to_physical_t3 (std::span<const double> tri, std::span<const double> ref)
{
  const std::size_t n = ref.size () / 2;
  std::vector<double> phy (2 * n);

  // Affine map from barycentric (1-r0-r1, r0, r1) to the physical triangle.
  for (std::size_t j = 0; j < n; ++j) {
    const double r0 = ref[2 * j];
    const double r1 = ref[2 * j + 1];
    for (int i = 0; i < 2; ++i)
      phy[2 * j + i] = tri[i] * (1.0 - r0 - r1) + tri[2 + i] * r0 + tri[4 + i] * r1;
  }

  return phy;
}

}  // namespace t8_mra

#endif
