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

}  // namespace t8_mra

#endif
