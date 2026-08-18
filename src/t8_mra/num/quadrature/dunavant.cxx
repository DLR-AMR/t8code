#ifdef T8_ENABLE_MRA

#include <cstddef>
#include <stdexcept>

#include "t8_mra/num/quadrature/dunavant.hxx"
#include "t8_mra/num/quadrature/dunavant_table.hxx"

namespace t8_mra
{

dunavant_quadrature
dunavant_rule (int rule)
{
  if (rule < 1 || rule > 20)
    throw std::out_of_range ("dunavant_rule: rule must be in [1, 20]");

  dunavant_quadrature q;

  const auto add = [&q] (double x, double y, double w) {
    q.points.push_back (x);
    q.points.push_back (y);
    q.weights.push_back (w);
  };

  // Each orbit expands into mult points by cyclic permutation of its
  // barycentric coordinates (x = b[k], y = b[k+1], third coord implied).
  for (const auto &o : dunavant_table::rule (rule)) {
    const auto &b = o.bary;
    switch (o.mult) {
    case 1:
      add (b[0], b[1], o.weight);
      break;
    case 3:
      for (int k = 0; k < 3; ++k)
        add (b[k % 3], b[(k + 1) % 3], o.weight);
      break;
    case 6:
      for (int k = 0; k < 3; ++k)
        add (b[k % 3], b[(k + 1) % 3], o.weight);
      for (int k = 0; k < 3; ++k)
        add (b[(k + 1) % 3], b[k % 3], o.weight);
      break;
    default:
      throw std::logic_error ("dunavant_rule: invalid orbit multiplicity");
    }
  }

  return q;
}

std::vector<double>
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
