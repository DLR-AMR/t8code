#pragma once

#include <vector>
#include "t8_eclass.h"

#ifdef T8_ENABLE_MRA

namespace t8_mra
{

template <t8_eclass TShape>
struct lmi_binary
{
  static constexpr int PATH_BITS = 0;
  static constexpr int LEVEL_BITS = 0;
  static constexpr int BASECELL_BITS = 0;
};

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
  const auto basecell = lmi.index & ((1ULL << BASECELL_BITS) - 1);
  lmi.index >>= BASECELL_BITS;

  // Extract level and remove level from lmi. Reduce refinement level by 1
  const auto level = (lmi.index & ((1u << LEVEL_BITS) - 1)) - 1;
  lmi.index >>= LEVEL_BITS;

  // Remove last path segment from lmi
  lmi.index >>= PATH_BITS;

  // Re-encode lmi with same basecell, decreased level and shorten path
  lmi.index = (lmi.index << (BASECELL_BITS + LEVEL_BITS)) | (level << BASECELL_BITS) | basecell;

  return lmi;
}

}  // namespace t8_mra

#endif
