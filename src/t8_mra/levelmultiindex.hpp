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

  levelmultiindex () = default;
  levelmultiindex (unsigned int level, size_t index);

  [[nodiscard]] unsigned int
  level () const noexcept;

  [[nodiscard]] size_t
  multiindex () const noexcept;

  [[nodiscard]] levelmultiindex
  parent (const levelmultiindex<TShape>& lmi) const;

  [[nodiscard]] std::vector<levelmultiindex>
  children (const levelmultiindex<TShape>& lmi) const noexcept;

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
[[nodiscard]] inline unsigned int
levelmultiindex<T8_ECLASS_TRIANGLE>::level () const noexcept
{
  return static_cast<unsigned int> ((index >> BASECELL_BITS >> LEVEL_BITS) & ((1ULL << LEVEL_BITS) - 1));
}

}  // namespace t8_mra

#endif
