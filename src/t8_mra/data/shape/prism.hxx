#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <cstddef>

#include "t8_eclass/t8_eclass.h"
#include "t8_element/t8_element.h"
#include "t8_schemes/t8_default/t8_default_prism/t8_dprism.h"
#include "t8_schemes/t8_default/t8_default_tri/t8_dtri.h"
#include "t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h"
#include "t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h"
#include "t8_schemes/t8_scheme.hxx"

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/data/levelmultiindex.hxx"
#include "t8_mra/data/triangle_order.hxx"

namespace t8_mra
{

template <>
struct lmi_properties<T8_ECLASS_PRISM>
{
  static constexpr int PATH_BITS = 3;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 20;
  static constexpr int NUM_CHILDREN = 8;
};

/// A prism is a triangle crossed with a line (t8_dprism), and refinement factors
/// the same way: t8code's child and ancestor ids are tri_id + 4 * line_id. Only
/// the triangle half carries a vertex order, so the Bey tracking below is the
/// triangle's, driven by the decomposed id.

/// Triangle child count, the radix of the prism child id.
inline constexpr int prism_tri_children = shape_traits<T8_ECLASS_TRIANGLE>::NUM_CHILDREN;

template <>
inline levelmultiindex<T8_ECLASS_PRISM>::levelmultiindex (size_t basecell, const t8_element_t *elem,
                                                          const t8_scheme *scheme) noexcept
  : index (basecell)
{
  std::array<int, 3> order = { 0, 1, 2 };

  const auto level = scheme->element_get_level (ECLASS, elem);
  const t8_dtri_t *tri = &reinterpret_cast<const t8_dprism_t *> (elem)->tri;
  t8_dtri_t ancestor;

  for (auto l = 0; l < level; ++l) {
    auto tmp = order;

    const auto prism_ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    const auto tri_ancestor_id = prism_ancestor_id % prism_tri_children;
    const auto line_ancestor_id = prism_ancestor_id / prism_tri_children;

    t8_dtri_ancestor (tri, l, &ancestor);
    triangle_order::invert_order (tmp);
    const auto tri_child_id = triangle_order::get_reference_children_order (ancestor.type, tri_ancestor_id, tmp);

    *this = jth_child (*this, tri_child_id + prism_tri_children * line_ancestor_id);
    triangle_order::get_point_order (order, t8_dtri_type_cid_to_beyid[ancestor.type][tri_ancestor_id]);
  }
}

template <>
inline std::array<int, 3>
levelmultiindex<T8_ECLASS_PRISM>::point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept
{
  std::array<int, 3> res = { 0, 1, 2 };
  const auto level = scheme->element_get_level (ECLASS, elem);
  const t8_dtri_t *tri = &reinterpret_cast<const t8_dprism_t *> (elem)->tri;
  t8_dtri_t ancestor;

  for (auto l = 0; l < level; ++l) {
    const auto tri_ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1) % prism_tri_children;

    t8_dtri_ancestor (tri, l, &ancestor);
    triangle_order::get_point_order (res, t8_dtri_type_cid_to_beyid[ancestor.type][tri_ancestor_id]);
  }

  return res;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
