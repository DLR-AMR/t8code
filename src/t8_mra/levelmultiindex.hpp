#pragma once

#include <vector>
#include "t8_eclass.h"

#ifdef T8_ENABLE_MRA

namespace t8_mra
{

/**
 * @brief Binary representation of each levelmultiindex. Has to be specialized
 * for each t8_eclass
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

  levelmultiindex ()
  {
    SC_ABORTF ("levelmultiindex has not been implemented for shape %d", TShape);
  };
  levelmultiindex (unsigned int level, size_t index);
  levelmultiindex (size_t index);

  [[nodiscard]] unsigned int
  level () const noexcept;

  [[nodiscard]] size_t
  multiindex () const noexcept;

  [[nodiscard]] static levelmultiindex
  parent (levelmultiindex<TShape> lmi);

  [[nodiscard]] static std::vector<levelmultiindex>
  children (levelmultiindex<TShape> lmi) noexcept;

 private:
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
levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex ()
{
  /// TODO
}

template <>
levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex (unsigned int level, size_t index)
{
  /// TODO
}

template <>
levelmultiindex<T8_ECLASS_TRIANGLE>::levelmultiindex (size_t _index): index (_index)
{
}

template <>
[[nodiscard]] inline unsigned int
levelmultiindex<T8_ECLASS_TRIANGLE>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS >> LEVEL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

/**
 * @brief Get parent of a given levelmultiindex.
 *
 * @param lmi Current levelmultiindex
 * @return Parent levelmultiindex
 */
template <>
[[nodiscard]] inline levelmultiindex<T8_ECLASS_TRIANGLE>
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
[[nodiscard]] inline std::vector<levelmultiindex<T8_ECLASS_TRIANGLE>>
levelmultiindex<T8_ECLASS_TRIANGLE>::children (levelmultiindex<T8_ECLASS_TRIANGLE> lmi) noexcept
{
  const auto NUM_CHILDREN = 4;

  std::vector<levelmultiindex<T8_ECLASS_TRIANGLE>> child_vec;
  child_vec.reserve (NUM_CHILDREN);

  // Gives path for the jth-child
  auto jth_path = [&] (size_t path, size_t j) -> size_t { return (path << PATH_BITS) | j; };

  // Extract basecell and remove basecell from lmi
  const size_t basecell = lmi.index & ((1u << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove basecell from lmi
  const size_t level = lmi.index & ((1u << LEVEL_BITS) - 1);
  lmi.index >>= LEVEL_BITS;

  // Construct all children: Same basecell, increase level by one, concat new
  // childpath
  for (size_t j = 0u; j < NUM_CHILDREN; ++j)
    child_vec.emplace_back ((jth_path (lmi.index, j) << (LEVEL_BITS + BASECELL_BITS)) | ((level + 1) << BASECELL_BITS)
                            | basecell);

  return child_vec;
}

// F R E E - F U N C T I O N S
template <typename TLmi>
inline TLmi
parent_lmi (TLmi lmi)
{
  return TLmi::parent (lmi);
}

template <typename TLmi>
inline std::vector<TLmi>
children_lmi (TLmi lmi)
{
  return TLmi::children (lmi);
}

}  // namespace t8_mra

#endif
