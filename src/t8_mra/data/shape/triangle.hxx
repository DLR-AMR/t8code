#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/triangle_order.hxx"

#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

namespace t8_mra
{

template <>
struct lmi_properties<T8_ECLASS_TRIANGLE>
{
  static constexpr int PATH_BITS = 2;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 21;
  static constexpr int NUM_CHILDREN = 4;
};

// The triangle vertex order changes with the Bey refinement type, so the element
// constructor and point_order_at_level track it down the ancestor chain.

template <>
inline levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex (size_t basecell, const t8_element_t *elem,
                                                             const t8_scheme *scheme) noexcept: index (basecell)
{
  std::array<int, 3> order = { 0, 1, 2 };
  const auto level = scheme->element_get_level (ECLASS, elem);
  t8_dtri_t ancestor;

  for (auto l = 0u; l < level; ++l) {
    auto tmp = order;

    const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    t8_dtri_ancestor ((t8_dtri_t *) elem, l, &ancestor);
    triangle_order::invert_order (tmp);
    const auto child_id = triangle_order::get_reference_children_order (ancestor.type, ancestor_id, tmp);

    *this = jth_child (*this, child_id);
    triangle_order::get_point_order (order, t8_dtri_type_cid_to_beyid[ancestor.type][ancestor_id]);
  }
}

template <>
inline std::array<int, 3>
levelmultiindex<T8_ECLASS_TRIANGLE>::point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept
{
  std::array<int, 3> res = { 0, 1, 2 };
  const auto level = scheme->element_get_level (ECLASS, elem);
  t8_dtri_t ancestor;

  for (auto l = 0; l < level; ++l) {
    const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    t8_dtri_ancestor ((t8_dtri_t *) elem, l, &ancestor);
    triangle_order::get_point_order (res, t8_dtri_type_cid_to_beyid[ancestor.type][ancestor_id]);
  }

  return res;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
