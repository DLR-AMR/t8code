#pragma once

#include <t8_mra/data/triangle_order.hpp>

#include <t8_eclass.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

#ifdef T8_ENABLE_MRA

namespace t8_mra
{

/**
 * @brief Properties of a levelmultiindex. Has to be specialized
 * for each t8_eclass.
 */
template <t8_eclass TShape>
struct lmi_properties
{
  static constexpr int PATH_BITS = 0;
  static constexpr int LEVEL_BITS = 0;
  static constexpr int BASECELL_BITS = 0;

  static constexpr int NUM_CHILDREN = 0;
};

/**
 * @brief Basic representation of levelmultiindex. Describes each cell in a grid
 * by its refinement level and path starting from a base cell.
 * Has to be specialized for each t8_eclass.
 * See: http://www.esaim-proc.org/10.1051/proc/201134003
 */
template <t8_eclass TShape>
struct levelmultiindex: public lmi_properties<TShape>
{

  static constexpr auto ECLASS = TShape;

  levelmultiindex () {
    // SC_ABORTF ("levelmultiindex has not been implemented for shape %d", TShape);
  };

  levelmultiindex (size_t _basecell) noexcept;
  levelmultiindex (size_t _basecell, const t8_element_t *elem, const t8_scheme *scheme) noexcept;

  bool
  operator== (const levelmultiindex &other) const noexcept
  {
    return index == other.index;
  }

  [[nodiscard]] unsigned int
  level () const noexcept;

  [[nodiscard]] size_t
  multiindex () const noexcept;

  /**
 * @brief Get parent of a given levelmultiindex.
 *
 * @param lmi Current levelmultiindex
 * @return Parent levelmultiindex
 */
  [[nodiscard]] static levelmultiindex
  parent (levelmultiindex<TShape> lmi);

  [[nodiscard]] static levelmultiindex
  jth_child (levelmultiindex<TShape> lmi, size_t j) noexcept;

  [[nodiscard]] static std::array<levelmultiindex, lmi_properties<TShape>::NUM_CHILDREN>
  children (levelmultiindex<TShape> lmi) noexcept;

  static std::array<int, 3>
  point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept;

  // private:
  size_t index;
};

// Specialization for TRIANGLE
template <>
struct lmi_properties<T8_ECLASS_TRIANGLE>
{
  static constexpr int PATH_BITS = 2;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 21;

  static constexpr int NUM_CHILDREN = 4;
};

// Specialization for LINE
template <>
struct lmi_properties<T8_ECLASS_LINE>
{
  static constexpr int PATH_BITS = 1;  // 2 children per refinement
  static constexpr int LEVEL_BITS = 6;
  static constexpr int BASECELL_BITS = 25;

  static constexpr int NUM_CHILDREN = 2;
};

// Specialization for QUAD
template <>
struct lmi_properties<T8_ECLASS_QUAD>
{
  static constexpr int PATH_BITS = 2;  // 4 children per refinement
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 21;

  static constexpr int NUM_CHILDREN = 4;
};

// Specialization for HEX
template <>
struct lmi_properties<T8_ECLASS_HEX>
{
  static constexpr int PATH_BITS = 3;  // 8 children per refinement
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 20;

  static constexpr int NUM_CHILDREN = 8;
};

// ============================================================================
// TRIANGLE SPECIALIZATIONS
// ============================================================================

template <>
inline levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex (size_t _basecell) noexcept: index (0u)
{
  index = (index << (LEVEL_BITS + BASECELL_BITS)) | _basecell;
}

template <>
inline levelmultiindex<T8_ECLASS_TRIANGLE>
levelmultiindex<T8_ECLASS_TRIANGLE>::jth_child (levelmultiindex<T8_ECLASS_TRIANGLE> lmi, size_t j) noexcept
{
  // Extract basecell and remove basecell from lmi
  const size_t basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove basecell from lmi
  const size_t level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  // Gives path for the jth-child
  const auto jth_path = (lmi.index << PATH_BITS) | j;

  // Construct all children: Same basecell, increase level by one, concat new
  // childpath
  lmi.index = (jth_path << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex (size_t _basecell, const t8_element_t *elem,
                                                             const t8_scheme *scheme) noexcept
  : levelmultiindex<T8_ECLASS_TRIANGLE> (_basecell)
{
  std::array<int, 3> order = { 0, 1, 2 };
  const auto level = scheme->element_get_level (ECLASS, elem);
  t8_dtri_t ancestor;

  for (auto l = 0u; l < level; ++l) {
    auto tmp = order;

    /// This should be possible in one function call...
    const auto ancestor_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    t8_dtri_ancestor ((t8_dtri_t *) elem, l, &ancestor);
    triangle_order::invert_order (tmp);
    const auto child_id = triangle_order::get_reference_children_order (ancestor.type, ancestor_id, tmp);

    *this = jth_child (*this, child_id);
    triangle_order::get_point_order (order, t8_dtri_type_cid_to_beyid[ancestor.type][ancestor_id]);
  }
}

template <>
inline unsigned int
levelmultiindex<T8_ECLASS_TRIANGLE>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

template <>
inline levelmultiindex<T8_ECLASS_TRIANGLE>
levelmultiindex<T8_ECLASS_TRIANGLE>::parent (levelmultiindex<T8_ECLASS_TRIANGLE> lmi)
{
#if T8_ENABLE_DEBUG
  if (parent.level == 0)
    SC_ABORTF ("levelmultiindices on level 0 do not have a parent %d", lmi.index);
#endif

  // Extract basecell and remove basecell from lmi
  const auto basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove level from lmi. Reduce refinement level by 1
  const auto level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  // Remove last path segment from lmi
  lmi.index >>= PATH_BITS;

  // Re-encode lmi with same basecell, decreased level by one and shorten path
  lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | ((level - 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline std::array<levelmultiindex<T8_ECLASS_TRIANGLE>, levelmultiindex<T8_ECLASS_TRIANGLE>::NUM_CHILDREN>
levelmultiindex<T8_ECLASS_TRIANGLE>::children (levelmultiindex<T8_ECLASS_TRIANGLE> lmi) noexcept
{
  std::array<levelmultiindex<T8_ECLASS_TRIANGLE>, levelmultiindex<T8_ECLASS_TRIANGLE>::NUM_CHILDREN> child_vec;

  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec[j] = jth_child (lmi, j);
  // child_vec[j] = jth_child (lmi.index, j);

  return child_vec;
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

// ============================================================================
// QUAD SPECIALIZATIONS
// ============================================================================

template <>
inline levelmultiindex<T8_ECLASS_QUAD>::levelmultiindex (size_t _basecell) noexcept: index (0u)
{
  index = (index << (LEVEL_BITS + BASECELL_BITS)) | _basecell;
}

template <>
inline levelmultiindex<T8_ECLASS_QUAD>
levelmultiindex<T8_ECLASS_QUAD>::jth_child (levelmultiindex<T8_ECLASS_QUAD> lmi, size_t j) noexcept
{
  // Extract basecell and remove basecell from lmi
  const size_t basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove level from lmi
  const size_t level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  // Append child ID to path
  const auto jth_path = (lmi.index << PATH_BITS) | j;

  // Reconstruct: Same basecell, increase level by one, append child path
  lmi.index = (jth_path << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline levelmultiindex<T8_ECLASS_QUAD>::levelmultiindex (size_t _basecell, const t8_element_t *elem,
                                                         const t8_scheme *scheme) noexcept
  : levelmultiindex<T8_ECLASS_QUAD> (_basecell)
{
  const auto level = scheme->element_get_level (ECLASS, elem);

  for (auto l = 0u; l < level; ++l) {
    const auto child_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    *this = jth_child (*this, child_id);
  }
}

template <>
inline unsigned int
levelmultiindex<T8_ECLASS_QUAD>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

template <>
inline levelmultiindex<T8_ECLASS_QUAD>
levelmultiindex<T8_ECLASS_QUAD>::parent (levelmultiindex<T8_ECLASS_QUAD> lmi)
{
  // Extract basecell and remove basecell from lmi
  const auto basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove level from lmi
  const auto level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  // Remove last path segment from lmi
  lmi.index >>= PATH_BITS;

  // Re-encode lmi with same basecell, decreased level by one and shortened path
  lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | ((level - 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline std::array<levelmultiindex<T8_ECLASS_QUAD>, levelmultiindex<T8_ECLASS_QUAD>::NUM_CHILDREN>
levelmultiindex<T8_ECLASS_QUAD>::children (levelmultiindex<T8_ECLASS_QUAD> lmi) noexcept
{
  std::array<levelmultiindex<T8_ECLASS_QUAD>, levelmultiindex<T8_ECLASS_QUAD>::NUM_CHILDREN> child_vec;

  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec[j] = jth_child (lmi, j);

  return child_vec;
}

template <>
inline std::array<int, 3>
levelmultiindex<T8_ECLASS_QUAD>::point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept
{
  // Cartesian elements don't need vertex ordering
  return { 0, 1, 2 };
}

// ============================================================================
// LINE SPECIALIZATIONS
// ============================================================================

template <>
inline levelmultiindex<T8_ECLASS_LINE>::levelmultiindex (size_t _basecell) noexcept: index (0u)
{
  index = (index << (LEVEL_BITS + BASECELL_BITS)) | _basecell;
}

template <>
inline levelmultiindex<T8_ECLASS_LINE>
levelmultiindex<T8_ECLASS_LINE>::jth_child (levelmultiindex<T8_ECLASS_LINE> lmi, size_t j) noexcept
{
  const size_t basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  const size_t level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  const auto jth_path = (lmi.index << PATH_BITS) | j;

  lmi.index = (jth_path << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline levelmultiindex<T8_ECLASS_LINE>::levelmultiindex (size_t _basecell, const t8_element_t *elem,
                                                         const t8_scheme *scheme) noexcept
  : levelmultiindex<T8_ECLASS_LINE> (_basecell)
{
  const auto level = scheme->element_get_level (ECLASS, elem);

  for (auto l = 0u; l < level; ++l) {
    const auto child_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    *this = jth_child (*this, child_id);
  }
}

template <>
inline unsigned int
levelmultiindex<T8_ECLASS_LINE>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

template <>
inline levelmultiindex<T8_ECLASS_LINE>
levelmultiindex<T8_ECLASS_LINE>::parent (levelmultiindex<T8_ECLASS_LINE> lmi)
{
  const auto basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  const auto level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  lmi.index >>= PATH_BITS;

  lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | ((level - 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline std::array<levelmultiindex<T8_ECLASS_LINE>, levelmultiindex<T8_ECLASS_LINE>::NUM_CHILDREN>
levelmultiindex<T8_ECLASS_LINE>::children (levelmultiindex<T8_ECLASS_LINE> lmi) noexcept
{
  std::array<levelmultiindex<T8_ECLASS_LINE>, levelmultiindex<T8_ECLASS_LINE>::NUM_CHILDREN> child_vec;

  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec[j] = jth_child (lmi, j);

  return child_vec;
}

template <>
inline std::array<int, 3>
levelmultiindex<T8_ECLASS_LINE>::point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept
{
  return { 0, 1, 2 };
}

// ============================================================================
// HEX SPECIALIZATIONS
// ============================================================================

template <>
inline levelmultiindex<T8_ECLASS_HEX>::levelmultiindex (size_t _basecell) noexcept: index (0u)
{
  index = (index << (LEVEL_BITS + BASECELL_BITS)) | _basecell;
}

template <>
inline levelmultiindex<T8_ECLASS_HEX>
levelmultiindex<T8_ECLASS_HEX>::jth_child (levelmultiindex<T8_ECLASS_HEX> lmi, size_t j) noexcept
{
  const size_t basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  const size_t level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  const auto jth_path = (lmi.index << PATH_BITS) | j;

  lmi.index = (jth_path << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline levelmultiindex<T8_ECLASS_HEX>::levelmultiindex (size_t _basecell, const t8_element_t *elem,
                                                        const t8_scheme *scheme) noexcept
  : levelmultiindex<T8_ECLASS_HEX> (_basecell)
{
  const auto level = scheme->element_get_level (ECLASS, elem);

  for (auto l = 0u; l < level; ++l) {
    const auto child_id = scheme->element_get_ancestor_id (ECLASS, elem, l + 1);
    *this = jth_child (*this, child_id);
  }
}

template <>
inline unsigned int
levelmultiindex<T8_ECLASS_HEX>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

template <>
inline levelmultiindex<T8_ECLASS_HEX>
levelmultiindex<T8_ECLASS_HEX>::parent (levelmultiindex<T8_ECLASS_HEX> lmi)
{
  const auto basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  const auto level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  lmi.index >>= PATH_BITS;

  lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | ((level - 1) << BASECELL_BITS) | basecell;

  return lmi;
}

template <>
inline std::array<levelmultiindex<T8_ECLASS_HEX>, levelmultiindex<T8_ECLASS_HEX>::NUM_CHILDREN>
levelmultiindex<T8_ECLASS_HEX>::children (levelmultiindex<T8_ECLASS_HEX> lmi) noexcept
{
  std::array<levelmultiindex<T8_ECLASS_HEX>, levelmultiindex<T8_ECLASS_HEX>::NUM_CHILDREN> child_vec;

  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec[j] = jth_child (lmi, j);

  return child_vec;
}

template <>
inline std::array<int, 3>
levelmultiindex<T8_ECLASS_HEX>::point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept
{
  return { 0, 1, 2 };
}

/// Levelmultiindex concept
template <typename T>
concept lmi_type = std::is_same_v<T, t8_mra::levelmultiindex<T::ECLASS>>;

// F R E E - F U N C T I O N S
template <lmi_type TLmi>
[[nodiscard]] inline TLmi
parent_lmi (TLmi lmi)
{
  return TLmi::parent (lmi);
}

template <lmi_type TLmi>
[[nodiscard]] inline TLmi
jth_child_lmi (TLmi lmi, size_t j)
{
  return TLmi::jth_child (lmi, j);
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
  using is_transparent = void;  // enable heterogeneous overloadsc
  using is_avalanching = void;

  size_t
  operator() (const TLmi &lmi) const
  {
    return lmi.index;
  }
};

}  // namespace std

#endif
