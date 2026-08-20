#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <array>

#include "t8_eclass/t8_eclass.h"
#include <t8_element/t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

namespace t8_mra
{

struct triangle_order
{
  static constexpr t8_eclass ECLASS = T8_ECLASS_TRIANGLE;

  using perm = std::array<int, 3>;

  static void
  get_point_order (perm &order, int cube_id)
  {
    const auto idx = permutation_index (order, point_perms);
    if (idx != -1)
      order = point_lookup[idx][cube_id];
  }

  static void
  invert_order (perm &order)
  {
    const auto idx = permutation_index (order, invert_perms);
    if (idx != -1)
      order = inverse_lookup[idx];
  }

  static void
  get_parent_order (perm &order)
  {
    order = parent_lookup[row_of (order)];
  }

  static int
  get_reference_children_order (int type, int child_id, const perm &order)
  {
    const auto &table = (type == 1) ? child_lookup_type_1 : child_lookup_type_2;
    return table[row_of (order)][child_id];
  }

  static void
  get_point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme, perm &order)
  {
    order = { 0, 1, 2 };
    const auto elem_level = scheme->element_get_level (ECLASS, elem);
    t8_dtri_t ancestor;

    for (auto l = 0; l < elem_level; ++l) {
      const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
      t8_dtri_ancestor ((t8_dtri_t *) elem, l, &ancestor);
      get_point_order (order, t8_dtri_type_cid_to_beyid[ancestor.type][ancestor_id]);
    }
  }

 private:
  /// Position of a permutation in the table, or -1 if absent.
  template <size_t N>
  static int
  permutation_index (const perm &order, const std::array<perm, N> &table)
  {
    const auto pos = std::ranges::find (table, order);
    return pos == table.end () ? -1 : static_cast<int> (pos - table.begin ());
  }

  /// Row index for the parent/children tables; unknown orders fall back to the last row.
  static int
  row_of (const perm &order)
  {
    const auto idx = permutation_index (order, point_perms);
    return idx < 0 ? 5 : idx;
  }

  static constexpr std::array<perm, 6> point_perms
    = { { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 } } };

  static constexpr std::array<perm, 6> invert_perms
    = { { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 2, 0 }, { 1, 0, 2 }, { 2, 0, 1 }, { 2, 1, 0 } } };

  static constexpr std::array<perm, 6> inverse_lookup
    = { { { 0, 1, 2 }, { 0, 2, 1 }, { 2, 0, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 1, 0 } } };

  static constexpr std::array<std::array<perm, 4>, 6> point_lookup = { {
    { { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 0, 2, 1 } } },
    { { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 2, 1, 0 } } },
    { { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 1, 0, 2 } } },
    { { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 2, 0, 1 } } },
    { { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 0, 1, 2 } } },
    { { { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 }, { 1, 2, 0 } } },
  } };

  static constexpr std::array<std::array<int, 4>, 6> child_lookup_type_1
    = { { { 1, 0, 2, 3 }, { 2, 0, 3, 1 }, { 3, 0, 1, 2 }, { 1, 0, 3, 2 }, { 2, 0, 1, 3 }, { 3, 0, 2, 1 } } };

  static constexpr std::array<std::array<int, 4>, 6> child_lookup_type_2
    = { { { 1, 2, 0, 3 }, { 2, 3, 0, 1 }, { 3, 1, 0, 2 }, { 1, 3, 0, 2 }, { 2, 1, 0, 3 }, { 3, 2, 0, 1 } } };

  static constexpr std::array<perm, 6> parent_lookup
    = { { { 1, 0, 2 }, { 0, 2, 1 }, { 2, 1, 0 }, { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } } };
};

}  // namespace t8_mra

#endif
