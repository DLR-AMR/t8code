#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include "t8_gtest_mra_forest.hxx"

namespace
{

using namespace mra_test;

constexpr double eps = 1e-12;

template <typename Cfg>
class mra_criteria: public ::testing::Test {};

TYPED_TEST_SUITE (mra_criteria, Configs, ConfigNames);


}  // namespace

#endif /* T8_ENABLE_MRA */
