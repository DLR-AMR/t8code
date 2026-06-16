#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_mra/data/levelmultiindex.hxx>
#include <array>
#include <type_traits>

namespace
{

template <t8_eclass TShape>
struct ShapeConfig
{
  static constexpr t8_eclass Shape = TShape;
};

using ShapeConfigs = ::testing::Types<ShapeConfig<T8_ECLASS_LINE>, ShapeConfig<T8_ECLASS_TRIANGLE>,
                                      ShapeConfig<T8_ECLASS_QUAD>, ShapeConfig<T8_ECLASS_HEX>>;

template <typename Config>
class mra_lmi: public ::testing::Test {};
TYPED_TEST_SUITE (mra_lmi, ShapeConfigs);

/* level() reads back the number of jth_child steps taken from a base cell. */
TYPED_TEST (mra_lmi, level_counts_refinement_steps)
{
  using lmi_t = t8_mra::levelmultiindex<TypeParam::Shape>;

  lmi_t lmi (7);  // arbitrary base cell, well inside BASECELL_BITS
  EXPECT_EQ (lmi.level (), 0u);

  for (auto l = 1u; l <= 5u; ++l) {
    lmi = lmi_t::jth_child (lmi, 0);
    EXPECT_EQ (lmi.level (), l);
  }
}

/* parent(jth_child(lmi, j)) == lmi for every child slot, and the base cell bits
 * survive the round trip. */
TYPED_TEST (mra_lmi, parent_inverts_jth_child)
{
  using lmi_t = t8_mra::levelmultiindex<TypeParam::Shape>;
  constexpr size_t basecell_mask = (static_cast<size_t> (1) << lmi_t::BASECELL_BITS) - 1;

  lmi_t lmi (13);
  for (auto step = 0; step < 4; ++step) {
    for (size_t j = 0; j < static_cast<size_t> (lmi_t::NUM_CHILDREN); ++j) {
      const auto child = lmi_t::jth_child (lmi, j);
      EXPECT_EQ (child.level (), lmi.level () + 1);
      EXPECT_EQ (child.index & basecell_mask, lmi.index & basecell_mask) << "base cell not preserved";
      EXPECT_EQ (lmi_t::parent (child).index, lmi.index) << "parent did not invert child " << j;
    }
    lmi = lmi_t::jth_child (lmi, 1);
  }
}

/* children() yields NUM_CHILDREN distinct cells, each one level finer and each
 * with lmi as parent; the free helpers agree with the static members. */
TYPED_TEST (mra_lmi, children_are_distinct_and_consistent)
{
  using lmi_t = t8_mra::levelmultiindex<TypeParam::Shape>;

  const lmi_t parent (5);
  const auto kids = lmi_t::children (parent);
  EXPECT_EQ (kids.size (), static_cast<size_t> (lmi_t::NUM_CHILDREN));

  for (size_t i = 0; i < kids.size (); ++i) {
    EXPECT_EQ (kids[i].level (), parent.level () + 1);
    EXPECT_EQ (lmi_t::parent (kids[i]).index, parent.index);
    for (size_t j = i + 1; j < kids.size (); ++j)
      EXPECT_NE (kids[i].index, kids[j].index) << "children " << i << " and " << j << " collide";
  }

  const auto free_kids = t8_mra::children_lmi (parent);
  for (size_t i = 0; i < kids.size (); ++i)
    EXPECT_EQ (free_kids[i].index, kids[i].index);
  EXPECT_EQ (t8_mra::parent_lmi (kids[0]).index, parent.index);
}

/* The lmi is its own hash (the value used as the dense-map key). */
TYPED_TEST (mra_lmi, hash_is_the_index)
{
  using lmi_t = t8_mra::levelmultiindex<TypeParam::Shape>;
  const auto lmi = lmi_t::jth_child (lmi_t (9), 1);
  EXPECT_EQ (std::hash<lmi_t> {}(lmi), lmi.index);
}

/* operator== compares the packed index (key equality in the dense map). */
TYPED_TEST (mra_lmi, equality_compares_the_index)
{
  using lmi_t = t8_mra::levelmultiindex<TypeParam::Shape>;
  const auto a = lmi_t::jth_child (lmi_t (6), 0);
  const auto same = lmi_t::jth_child (lmi_t (6), 0);
  const auto other = lmi_t::jth_child (lmi_t (6), 1);

  EXPECT_TRUE (a == same);
  EXPECT_FALSE (a == other);
}


}  // namespace

#endif  // T8_ENABLE_MRA
