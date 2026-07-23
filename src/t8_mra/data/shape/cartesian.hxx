#pragma once

#ifdef T8_ENABLE_MRA

#include "t8_mra/data/levelmultiindex.hxx"

namespace t8_mra
{

template <>
struct lmi_properties<T8_ECLASS_LINE>
{
  static constexpr int PATH_BITS = 1;
  static constexpr int LEVEL_BITS = 6;
  static constexpr int BASECELL_BITS = 25;
  static constexpr int NUM_CHILDREN = 2;
};

template <>
struct lmi_properties<T8_ECLASS_QUAD>
{
  static constexpr int PATH_BITS = 2;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 21;
  static constexpr int NUM_CHILDREN = 4;
};

template <>
struct lmi_properties<T8_ECLASS_HEX>
{
  static constexpr int PATH_BITS = 3;
  static constexpr int LEVEL_BITS = 5;
  static constexpr int BASECELL_BITS = 20;
  static constexpr int NUM_CHILDREN = 8;
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA
