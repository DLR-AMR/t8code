#pragma once

#ifdef T8_ENABLE_MRA

#include <array>

#include "t8_eclass/t8_eclass.h"

namespace t8_mra
{

namespace detail
{

using triangle_perm = std::array<int, 3>;

/// Perfect hash of a permutation of {0,1,2} into [0, 27).
[[nodiscard]] constexpr int
triangle_perm_code (const triangle_perm &order) noexcept
{
  return 9 * order[0] + 3 * order[1] + order[2];
}

/// Row of each permutation in table, indexed by triangle_perm_code; -1 where absent.
template <size_t N>
[[nodiscard]] constexpr std::array<signed char, 27>
triangle_perm_rows (const std::array<triangle_perm, N> &table) noexcept
{
  std::array<signed char, 27> rows;
  rows.fill (-1);

  for (auto row = 0u; row < N; ++row)
    rows[triangle_perm_code (table[row])] = static_cast<signed char> (row);

  return rows;
}

}  // namespace detail

struct triangle_order
{
  static constexpr t8_eclass ECLASS = T8_ECLASS_TRIANGLE;

  using perm = std::array<int, 3>;

  static void
  get_point_order (perm &order, int cube_id)
  {
    const auto row = point_rows[detail::triangle_perm_code (order)];
    if (row >= 0)
      order = point_lookup[row][cube_id];
  }

  static void
  invert_order (perm &order)
  {
    const auto row = invert_rows[detail::triangle_perm_code (order)];
    if (row >= 0)
      order = inverse_lookup[row];
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

 private:
  /// Row index for the parent/children tables; unknown orders fall back to the last row.
  static int
  row_of (const perm &order)
  {
    const auto row = point_rows[detail::triangle_perm_code (order)];
    return row < 0 ? 5 : row;
  }

  static constexpr std::array<perm, 6> point_perms
    = { { { 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }, { 0, 2, 1 }, { 1, 0, 2 }, { 2, 1, 0 } } };

  static constexpr std::array<perm, 6> invert_perms
    = { { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 2, 0 }, { 1, 0, 2 }, { 2, 0, 1 }, { 2, 1, 0 } } };

  static constexpr auto point_rows = detail::triangle_perm_rows (point_perms);
  static constexpr auto invert_rows = detail::triangle_perm_rows (invert_perms);

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
