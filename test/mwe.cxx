#include <gtest/gtest.h>

class class_mwe: public testing::TestWithParam<int> {
 protected:
  void
  SetUp () override
  {
  }
  void
  TearDown () override
  {
  }
};

TEST_P (class_mwe, mwe)
{
  return;
}

INSTANTIATE_TEST_SUITE_P (mwe, class_mwe, testing::Range (0, 1));