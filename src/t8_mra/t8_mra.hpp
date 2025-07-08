#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_eclass.h"

#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/data/levelindex_map.hpp"

namespace t8_mra
{

template <t8_eclass TShape>
class multiscale_types {
  static constexpr unsigned short DIM = 0;
};

template <>
class multiscale_types<T8_ECLASS_TRIANGLE> {
  static constexpr unsigned short DIM = 2;
};

template <t8_eclass TShape, int U, int ORDER>
class multiscale: multiscale_types<TShape> {
  using element_t = data_per_element<TShape, U, ORDER>;
  using levelmultiindex = levelmultiindex<TShape>;

  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = element_t::P_DIM;

  void
  multiscale_transformation (t8_mra::levelindex_map<element_t>& grid_hierarchy, unsigned int l_min, unsigned int l_max)
  {
    for (auto l = l_max; l > l_min; --l) {
      for (const auto& [lmi, val] : grid_hierarchy.level_map[l]) {
        const auto parent_lmi = t8_mra::parent_lmi (lmi);

        if (grid_hierarchy.contains (l - 1, parent_lmi))
          continue;

        const auto children = t8_mra::children_lmi (parent_lmi);
        std::array<element_t, levelmultiindex::NUM_CHILDREN> child_data;
        element_t parent_data;

        for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k)
          child_data[k] = grid_hierarchy.get (l, children[k]);

        parent_data.order = child_data[0].order;
        triangle_order::get_parent_order (parent_data.order);

        for (auto i = 0u; i < P_DIM; ++i) {
          auto u_sum = 0.0;
          auto d_sum = 0.0;

          for (auto j = 0u; j < P_DIM; ++j) {
            for (auto k = 0u; k < levelmultiindex::NUM_CHILDREN; ++k) {
              /// TODO mask coeffs
            }
          }
        }
      }
    }
  }
};

}  // namespace t8_mra

#endif
