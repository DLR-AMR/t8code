#pragma once

#include <t8_eclass/t8_eclass.h>
#include <t8_element/t8_element.h>
#include <t8_schemes/t8_scheme.hxx>

#include <array>

#ifdef T8_ENABLE_MRA

namespace t8_mra
{

// ============================================================================
// Adding a new element shape
// ============================================================================
// In data/shape/<shape>.hxx:
// 1. Specialize lmi_properties<Shape> with the bit widths (PATH_BITS =
//    log2(NUM_CHILDREN), LEVEL_BITS, BASECELL_BITS; their sum must be <= 64)
//    and NUM_CHILDREN. That alone gives a working lmi for any shape whose
//    children carry no vertex-ordering information (all cartesian shapes).
// 2. Only if the shape's vertex order changes under refinement (like the
//    simplex Bey refinement): specialize the element constructor and
//    point_order_at_level (see data/shape/triangle.hxx).

/// Bit layout (path | level | basecell) and child count of an lmi, per shape.
template <t8_eclass TShape>
struct lmi_properties
{
  static constexpr int PATH_BITS = 0;
  static constexpr int LEVEL_BITS = 0;
  static constexpr int BASECELL_BITS = 0;
  static constexpr int NUM_CHILDREN = 0;
};

/**
 * @brief Cell identifier: (base tree, level, path) packed into one size_t.
 *
 * Describes any cell of an adaptive grid by its refinement level and the
 * child-id path from its base tree, independently of which rank holds it.
 * The bit layout and child count come from lmi_properties<TShape>; the
 * index arithmetic below is shape-independent. Only the construction from a
 * forest element and the vertex order can be shape-specific (see TRIANGLE).
 * See http://www.esaim-proc.org/10.1051/proc/201134003
 */
template <t8_eclass TShape>
struct levelmultiindex
{
  static constexpr auto ECLASS = TShape;
  static constexpr int PATH_BITS = lmi_properties<TShape>::PATH_BITS;
  static constexpr int LEVEL_BITS = lmi_properties<TShape>::LEVEL_BITS;
  static constexpr int BASECELL_BITS = lmi_properties<TShape>::BASECELL_BITS;
  static constexpr int NUM_CHILDREN = lmi_properties<TShape>::NUM_CHILDREN;

  size_t index;

  levelmultiindex () = default;

  levelmultiindex (size_t basecell) noexcept: index (basecell)
  {
  }

  /// Construct from a forest element by walking its ancestor child-ids down to
  /// its level. Generic (cartesian) form; TRIANGLE specializes this to also
  /// track the vertex order.
  levelmultiindex (size_t basecell, const t8_element_t *elem, const t8_scheme *scheme) noexcept: index (basecell)
  {
    const auto level = scheme->element_get_level (ECLASS, elem);
    for (auto l = 0u; l < level; ++l)
      *this = jth_child (*this, scheme->element_get_ancestor_id (ECLASS, elem, l + 1));
  }

  bool
  operator== (const levelmultiindex &other) const noexcept
  {
    return index == other.index;
  }

  [[nodiscard]] unsigned int
  level () const noexcept
  {
    return static_cast<unsigned int> ((index >> BASECELL_BITS) & ((1ULL << LEVEL_BITS) - 1));
  }

  /// Child j: append j to the path, increment the level.
  [[nodiscard]] static levelmultiindex
  jth_child (levelmultiindex lmi, size_t j) noexcept
  {
    const size_t basecell = lmi.index & ((1ULL << BASECELL_BITS) - 1);
    lmi.index >>= BASECELL_BITS;
    const size_t level = lmi.index & ((1ULL << LEVEL_BITS) - 1);
    lmi.index >>= LEVEL_BITS;

    const auto jth_path = (lmi.index << PATH_BITS) | j;
    lmi.index = (jth_path << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS) | basecell;

    return lmi;
  }

  /// Parent: drop the last path segment, decrement the level.
  [[nodiscard]] static levelmultiindex
  parent (levelmultiindex lmi)
  {
#if T8_ENABLE_DEBUG
    if (lmi.level () == 0)
      SC_ABORTF ("levelmultiindices on level 0 do not have a parent %zu", lmi.index);
#endif
    const auto basecell = lmi.index & ((1ULL << BASECELL_BITS) - 1);
    lmi.index >>= BASECELL_BITS;
    const auto level = lmi.index & ((1ULL << LEVEL_BITS) - 1);
    lmi.index >>= LEVEL_BITS;
    lmi.index >>= PATH_BITS;

    lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | ((level - 1) << BASECELL_BITS) | basecell;

    return lmi;
  }

  [[nodiscard]] static std::array<levelmultiindex, NUM_CHILDREN>
  children (levelmultiindex lmi) noexcept
  {
    std::array<levelmultiindex, NUM_CHILDREN> child_vec;
    for (size_t j = 0u; j < NUM_CHILDREN; ++j)
      child_vec[j] = jth_child (lmi, j);

    return child_vec;
  }

  /// Reference vertex order at the element's level. Cartesian elements need no
  /// reordering (identity); TRIANGLE specializes.
  static std::array<int, 3>
  point_order_at_level (const t8_element_t * /*elem*/, const t8_scheme * /*scheme*/) noexcept
  {
    return { 0, 1, 2 };
  }
};

/// Concept: an lmi is a levelmultiindex of its own ECLASS.
template <typename T>
concept lmi_type = std::is_same_v<T, t8_mra::levelmultiindex<T::ECLASS>>;

template <lmi_type TLmi>
[[nodiscard]] inline TLmi
parent_lmi (TLmi lmi)
{
  return TLmi::parent (lmi);
}

template <lmi_type TLmi>
[[nodiscard]] inline std::array<TLmi, TLmi::NUM_CHILDREN>
children_lmi (TLmi lmi)
{
  return TLmi::children (lmi);
}

}  // namespace t8_mra

namespace std
{

template <t8_mra::lmi_type TLmi>
struct hash<TLmi>
{
  using is_transparent = void;

  size_t
  operator() (const TLmi &lmi) const
  {
    return lmi.index;
  }
};

}  // namespace std

// Per-shape specializations (defined after the primary template).
#include "t8_mra/data/shape/cartesian.hxx"
#include "t8_mra/data/shape/triangle.hxx"

#endif
