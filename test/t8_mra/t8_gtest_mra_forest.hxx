#ifndef T8_GTEST_MRA_FOREST_HXX
#define T8_GTEST_MRA_FOREST_HXX

#ifdef T8_ENABLE_MRA

#include <gtest/gtest.h>

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_mra/t8_mra.hxx>

#include <array>
#include <cmath>
#include <string>

namespace mra_test
{

template <t8_eclass TShape, int U_, int P_>
struct Config
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int U = U_;
  static constexpr int P = P_;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;
};

using Configs
  = ::testing::Types<Config<T8_ECLASS_TRIANGLE, 1, 1>, Config<T8_ECLASS_TRIANGLE, 1, 2>,
                     Config<T8_ECLASS_TRIANGLE, 1, 3>, Config<T8_ECLASS_TRIANGLE, 1, 4>,
                     Config<T8_ECLASS_TRIANGLE, 2, 2>, Config<T8_ECLASS_QUAD, 1, 1>, Config<T8_ECLASS_QUAD, 1, 2>,
                     Config<T8_ECLASS_QUAD, 1, 3>, Config<T8_ECLASS_QUAD, 1, 4>, Config<T8_ECLASS_QUAD, 2, 2>,
                     Config<T8_ECLASS_HEX, 1, 2>, Config<T8_ECLASS_HEX, 1, 3>, Config<T8_ECLASS_HEX, 2, 2>>;

struct ConfigNames
{
  template <typename T>
  static std::string
  GetName (int)
  {
    return std::string (t8_eclass_to_string[T::Shape]) + "_U" + std::to_string (T::U) + "_P" + std::to_string (T::P);
  }
};

/* Smooth test function */

}  // namespace mra_test

#endif /* T8_ENABLE_MRA */

#endif /* T8_GTEST_MRA_FOREST_HXX */
