#pragma once

#include <vector>

#include <t8_mra/triangle_order.hpp>

#include <t8_eclass.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

#ifdef T8_ENABLE_MRA

namespace t8_mra
{

/**
 * @brief Binary representation of each levelmultiindex. Has to be specialized
 * for each t8_eclass.
 */
template <t8_eclass TShape>
struct lmi_binary
{
  static constexpr int PATH_BITS = 0;
  static constexpr int LEVEL_BITS = 0;
  static constexpr int BASECELL_BITS = 0;
};

/**
 * @brief Basic representation of levelmultiindex. Describes each cell in a grid
 * by its refinement level and path starting from a base cell.
 * Has to be specialized for each t8_eclass.
 * See: http://www.esaim-proc.org/10.1051/proc/201134003
 */
template <t8_eclass TShape>
struct levelmultiindex: public lmi_binary<TShape>
{

  static constexpr auto ECLASS = TShape;

  levelmultiindex ()
  {
    SC_ABORTF ("levelmultiindex has not been implemented for shape %d", TShape);
  };
  levelmultiindex (size_t _basecell) noexcept;
  levelmultiindex (size_t _basecell, const t8_element_t *elem, const t8_scheme *scheme) noexcept;

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

  [[nodiscard]] static std::vector<levelmultiindex>
  children (levelmultiindex<TShape> lmi) noexcept;

  static std::array<int, 3>
  point_order_at_level (const t8_element_t *elem, const t8_scheme *scheme) noexcept;

  // private:
  size_t index;
};

template <>
struct lmi_binary<T8_ECLASS_TRIANGLE>
{
  static constexpr int PATH_BITS = 2;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 21;
};

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
  return static_cast<unsigned int> ((index >> BASECELL_BITS >> LEVEL_BITS) & ((1ULL << LEVEL_BITS) - 1));
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
inline std::vector<levelmultiindex<T8_ECLASS_TRIANGLE>>
levelmultiindex<T8_ECLASS_TRIANGLE>::children (levelmultiindex<T8_ECLASS_TRIANGLE> lmi) noexcept
{
  const auto NUM_CHILDREN = 4;

  std::vector<levelmultiindex<T8_ECLASS_TRIANGLE>> child_vec;
  child_vec.reserve (NUM_CHILDREN);

  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec.emplace_back (jth_child (lmi.index, j));

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

// F R E E - F U N C T I O N S
template <typename TLmi>
[[nodiscard]] inline TLmi
parent_lmi (TLmi lmi)
{
  return TLmi::parent (lmi);
}

template <typename TLmi>
[[nodiscard]] inline TLmi
jth_child_lmi (TLmi lmi, size_t j)
{
  return TLmi::jth_child (lmi, j);
}

template <typename TLmi>
[[nodiscard]] inline std::vector<TLmi>
children_lmi (TLmi lmi)
{
  return TLmi::children (lmi);
}

}  // namespace t8_mra

#endif
