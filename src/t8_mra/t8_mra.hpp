#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_eclass.h"

#include "t8_mra/data/cell_data.hpp"
#include "t8_mra/data/levelmultiindex.hpp"
#include "t8_mra/data/levelindex_map.hpp"
#include "t8_mra/num/mask_coefficients.hpp"

namespace t8_mra
{

template <t8_eclass TShape>
struct multiscale_data
{
  static constexpr unsigned short DIM = 0;

  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

template <>
struct multiscale_data<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned short DIM = 2;

  std::vector<t8_mra::mat> mask_coefficients;
  std::vector<t8_mra::mat> inverse_mask_coefficients;
};

/// TODO naming P -> ORDER?
template <t8_eclass TShape, int U, int P>
class multiscale: public multiscale_data<TShape> {
  using element_t = data_per_element<TShape, U, P>;
  using levelmultiindex = levelmultiindex<TShape>;

  static constexpr unsigned int DIM = element_t::DIM;
  static constexpr unsigned int U_DIM = U;
  static constexpr unsigned int P_DIM = P;

  static constexpr unsigned int DOF = element_t::DOF;
  static constexpr unsigned int W_DOF = element_t::W_DOF;

 public:
  multiscale ()
  {
    t8_mra::initialize_mask_coefficients<TShape> (P_DIM, DOF, multiscale_data<TShape>::mask_coefficients,
                                                  multiscale_data<TShape>::inverse_mask_coefficients);
  }

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

        for (auto i = 0u; i < DOF; ++i) {
          auto u_sum = 0.0;
          auto d_sum = 0.0;

          for (auto j = 0u; j < DOF; ++j) {
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
