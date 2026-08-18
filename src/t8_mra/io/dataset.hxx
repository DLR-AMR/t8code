#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace t8_mra
{

/**
 * @brief Regular-grid sampler usable as an MRA projection func
 *
 * Holds U components on a DIM-dimensional structured grid and multilinearly
 * interpolates them. The call operators take DIM positional doubles and return
 * std::array<double, U>, matching the func signature initialize_data expects,
 * so a dataset drops straight into initialize_data / initialize_data_adaptive.
 * Points outside the grid box clamp to the boundary.
 *
 * Value layout: node-major with axis 0 slowest (axis DIM-1 fastest), the U
 * components contiguous per node.
 */
template <unsigned int DIM, unsigned int U>
class structured_dataset {
 public:
  structured_dataset () = default;

  structured_dataset (const std::array<std::size_t, DIM> &num_nodes, const std::array<double, DIM> &origin,
                      const std::array<double, DIM> &spacing, std::vector<double> values)
    : n (num_nodes), origin (origin), spacing (spacing), values (std::move (values))
  {
    std::size_t total = U;
    for (auto d = 0u; d < DIM; ++d) {
      if (n[d] == 0)
        throw std::invalid_argument ("structured_dataset: zero nodes on an axis");
      total *= n[d];
    }
    if (this->values.size () != total)
      throw std::invalid_argument ("structured_dataset: value count does not match grid size * U");

    stride[DIM - 1] = 1;
    for (auto d = DIM - 1; d-- > 0;)
      stride[d] = stride[d + 1] * n[d + 1];
  }

  [[nodiscard]] std::array<double, U>
  operator() (double x) const
    requires (DIM == 1)
  {
    return sample ({ x });
  }

  [[nodiscard]] std::array<double, U>
  operator() (double x, double y) const
    requires (DIM == 2)
  {
    return sample ({ x, y });
  }

  [[nodiscard]] std::array<double, U>
  operator() (double x, double y, double z) const
    requires (DIM == 3)
  {
    return sample ({ x, y, z });
  }

 private:
  /// Multilinear interpolation at p, clamped to the grid box.
  [[nodiscard]] std::array<double, U>
  sample (const std::array<double, DIM> &p) const
  {
    std::array<std::size_t, DIM> base;
    std::array<double, DIM> t;
    for (auto d = 0u; d < DIM; ++d) {
      if (n[d] == 1) {
        base[d] = 0;
        t[d] = 0.0;
        continue;
      }
      const auto local = std::clamp ((p[d] - origin[d]) / spacing[d], 0.0, static_cast<double> (n[d] - 1));
      const auto i = static_cast<std::size_t> (std::clamp (std::floor (local), 0.0, static_cast<double> (n[d] - 2)));
      base[d] = i;
      t[d] = local - static_cast<double> (i);
    }

    std::array<double, U> out = {};
    for (unsigned int corner = 0; corner < (1u << DIM); ++corner) {
      double w = 1.0;
      std::size_t flat = 0;
      for (auto d = 0u; d < DIM; ++d) {
        const bool upper = corner & (1u << d);
        flat += (base[d] + (upper ? 1 : 0)) * stride[d];
        w *= upper ? t[d] : 1.0 - t[d];
      }
      const auto *node = &values[flat * U];
      for (auto u = 0u; u < U; ++u)
        out[u] += w * node[u];
    }
    return out;
  }

  std::array<std::size_t, DIM> n = {};
  std::array<std::size_t, DIM> stride = {};
  std::array<double, DIM> origin = {};
  std::array<double, DIM> spacing = {};
  std::vector<double> values;
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
