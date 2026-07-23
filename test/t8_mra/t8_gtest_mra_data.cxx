#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/data/levelmultiindex.hxx>
#include <t8_mra/data/levelindex_map.hxx>
#include <t8_mra/data/levelindex_set.hxx>
#include <t8_mra/data/element_data.hxx>

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

/* ---- containers: shape-independent, one representative shape ---- */

using lmi_q = t8_mra::levelmultiindex<T8_ECLASS_QUAD>;

TEST (mra_levelindex_map, insert_get_contains_erase)
{
  t8_mra::levelindex_map<lmi_q, int> map (6);

  const auto a = lmi_q::jth_child (lmi_q (2), 0);  // level 1
  const auto b = lmi_q::jth_child (a, 3);          // level 2

  EXPECT_FALSE (map.contains (a));
  EXPECT_EQ (map.find (a), nullptr);

  map.insert (a, 11);
  map.insert (b, 22);

  EXPECT_TRUE (map.contains (a));
  EXPECT_EQ (map.get (a), 11);
  EXPECT_EQ (*map.find (b), 22);
  EXPECT_EQ (map.size (), 2u);
  EXPECT_EQ (map.size (a.level ()), 1u);
  EXPECT_EQ (map.size (b.level ()), 1u);

  /// in-place mutation through find
  *map.find (a) = 99;
  EXPECT_EQ (map.get (a), 99);

  /// overwrite via insert
  map.insert (a, 5);
  EXPECT_EQ (map.get (a), 5);
  EXPECT_EQ (map.size (), 2u);

  map.erase (a);
  EXPECT_FALSE (map.contains (a));
  EXPECT_EQ (map.size (), 1u);

  map.erase_all ();
  EXPECT_EQ (map.size (), 0u);
}

TEST (mra_levelindex_map, erase_level_clears_only_that_level)
{
  t8_mra::levelindex_map<lmi_q, int> map (6);
  const auto l1a = lmi_q::jth_child (lmi_q (0), 0);
  const auto l1b = lmi_q::jth_child (lmi_q (0), 1);
  const auto l2 = lmi_q::jth_child (l1a, 0);

  map.insert (l1a, 1);
  map.insert (l1b, 2);
  map.insert (l2, 3);

  map.erase (l1a.level ());
  EXPECT_EQ (map.size (l1a.level ()), 0u);
  EXPECT_EQ (map.size (l2.level ()), 1u);
  EXPECT_TRUE (map.contains (l2));
}

/* iterating through a refinement level */
TEST (mra_levelindex_map, level_view_iterates_that_level)
{
  t8_mra::levelindex_map<lmi_q, int> map (6);

  const auto root = lmi_q (3);
  const auto kids = lmi_q::children (root);

  for (size_t i = 0; i < kids.size (); ++i)
    map.insert (kids[i], static_cast<int> (10 + i));

  map.insert (lmi_q::jth_child (kids[0], 0), 99);  // a level-2 cell

  const unsigned int level = 1;
  EXPECT_EQ (map[level].size (), kids.size ());

  size_t counted = 0;
  for (auto it = map.begin (level); it != map.end (level); ++it) {
    EXPECT_EQ (it->first.level (), level);
    EXPECT_TRUE (map.contains (it->first));
    ++counted;
  }
  EXPECT_EQ (counted, kids.size ());
}

/* The (level, key) overloads agree with the lmi overloads. */
TEST (mra_levelindex_map, level_key_overloads_match_lmi_overloads)
{
  t8_mra::levelindex_map<lmi_q, int> map (6);
  const auto a = lmi_q::jth_child (lmi_q (4), 2);

  map.insert (a.level (), a.index, 7);
  EXPECT_TRUE (map.contains (a.level (), a.index));
  EXPECT_TRUE (map.contains (a));
  EXPECT_EQ (map.get (a.level (), a.index), 7);
  EXPECT_EQ (map.get (a), 7);

  map.erase (a.level (), a.index);
  EXPECT_FALSE (map.contains (a));
}

TEST (mra_levelindex_set, insert_contains_erase_sizes)
{
  t8_mra::levelindex_set<lmi_q> set (6);
  const auto a = lmi_q::jth_child (lmi_q (1), 2);
  const auto b = lmi_q::jth_child (a, 1);

  set.insert (a);
  set.insert (b);
  set.insert (a);  // idempotent

  EXPECT_TRUE (set.contains (a));
  EXPECT_TRUE (set.contains (b));
  EXPECT_EQ (set.size (), 2u);
  EXPECT_EQ (set.size (a.level ()), 1u);

  set.erase (b);
  EXPECT_FALSE (set.contains (b));
  EXPECT_EQ (set.size (), 1u);

  set.erase_all ();
  EXPECT_EQ (set.size (), 0u);
}

TEST (mra_levelindex_set, level_view_and_level_key_overloads)
{
  t8_mra::levelindex_set<lmi_q> set (6);
  const auto kids = lmi_q::children (lmi_q (2));
  for (const auto &k : kids)
    set.insert (k);

  const unsigned int level = 1;
  EXPECT_EQ (set[level].size (), kids.size ());

  size_t counted = 0;
  for (auto it = set.begin (level); it != set.end (level); ++it)
    ++counted;
  EXPECT_EQ (counted, kids.size ());

  const auto a = kids[0];
  EXPECT_TRUE (set.contains (a.level (), a.index));
  set.erase (a.level (), a.index);
  EXPECT_FALSE (set.contains (a));
}

/* ---- element_data ---- */

TEST (mra_element_data, dg_idx_is_a_bijection)
{
  using elem = t8_mra::element_data<T8_ECLASS_QUAD, 3, 4>;  // U=3, DOF = 4*4 = 16
  constexpr auto U = elem::U_DIM;
  constexpr auto DOF = elem::DOF;

  std::array<bool, U * DOF> seen {};
  for (size_t u = 0; u < U; ++u)
    for (size_t p = 0; p < DOF; ++p) {
      const auto idx = elem::dg_idx (u, p);
      ASSERT_LT (idx, U * DOF);
      EXPECT_FALSE (seen[idx]) << "dg_idx collision at (" << u << "," << p << ")";
      seen[idx] = true;
    }
}

TEST (mra_element_data, wavelet_idx_is_a_bijection)
{
  using det = t8_mra::detail_data<T8_ECLASS_TRIANGLE, 2, 3>;
  constexpr auto U = det::U_DIM;
  constexpr auto DOF = det::DOF;
  constexpr auto NC = det::NUM_CHILDREN;

  std::array<bool, NC * U * DOF> seen {};
  for (size_t k = 0; k < NC; ++k)
    for (size_t u = 0; u < U; ++u)
      for (size_t p = 0; p < DOF; ++p) {
        const auto idx = det::wavelet_idx (k, u, p);
        ASSERT_LT (idx, NC * U * DOF);
        EXPECT_FALSE (seen[idx]) << "wavelet_idx collision at (" << k << "," << u << "," << p << ")";
        seen[idx] = true;
      }
}

/* The repartition/ghost wire ships element_data as raw bytes. */
TEST (mra_element_data, is_trivially_copyable)
{
  EXPECT_TRUE ((std::is_trivially_copyable_v<t8_mra::element_data<T8_ECLASS_HEX, 2, 2>>) );
  EXPECT_TRUE ((std::is_trivially_copyable_v<t8_mra::detail_data<T8_ECLASS_QUAD, 1, 3>>) );
}

/* triangle_order: TRIANGLE-specific Bey vertex-order tables (pure) */

/* The six permutations of three vertices. */
const std::array<std::array<int, 3>, 6> all_orders
  = { { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } } };

template <std::size_t N>
bool
is_permutation (const std::array<int, N> &p)
{
  std::array<bool, N> seen {};
  for (int v : p) {
    if (v < 0 || v >= static_cast<int> (N) || seen[v])
      return false;
    seen[v] = true;
  }
  return true;
}

/* invert_order produces the true inverse permutation q (q[p[i]] = i). */
TEST (mra_triangle_order, invert_order_is_the_permutation_inverse)
{
  for (const auto &order : all_orders) {
    std::array<int, 3> expected {};
    for (int i = 0; i < 3; ++i)
      expected[order[i]] = i;

    auto tmp_order = order;
    t8_mra::triangle_order::invert_order (tmp_order);
    EXPECT_EQ (tmp_order, expected) << "wrong inverse for {" << order[0] << order[1] << order[2] << "}";
  }
}

/* Inverting twice is the identity. */
TEST (mra_triangle_order, invert_order_is_an_involution)
{
  for (const auto &order : all_orders) {
    auto tmp_order = order;
    t8_mra::triangle_order::invert_order (tmp_order);
    t8_mra::triangle_order::invert_order (tmp_order);
    EXPECT_EQ (tmp_order, order);
  }
}

/* get_point_order and get_parent_order map a permutation to a permutation; the
 * identity order with child 0 stays the identity. */
TEST (mra_triangle_order, order_maps_stay_permutations)
{
  for (const auto &order : all_orders)
    for (int bey = 0; bey < 4; ++bey) {
      auto tmp_order = order;
      t8_mra::triangle_order::get_point_order (tmp_order, bey);
      EXPECT_TRUE (is_permutation<3> (tmp_order)) << "get_point_order broke the permutation";
    }

  std::array<int, 3> id { 0, 1, 2 };
  t8_mra::triangle_order::get_point_order (id, 0);
  EXPECT_EQ ((std::array<int, 3> { 0, 1, 2 }), id);

  for (const auto &order : all_orders) {
    auto tmp_order = order;
    t8_mra::triangle_order::get_parent_order (tmp_order);
    EXPECT_TRUE (is_permutation<3> (tmp_order)) << "get_parent_order broke the permutation";
  }
}

/* For a fixed parent order and Bey type, child_id -> reference index is a
 * bijection over {0,1,2,3} (the four children land on four distinct slots). */
TEST (mra_triangle_order, reference_children_order_is_a_bijection)
{
  for (int type = 1; type <= 2; ++type)
    for (const auto &order : all_orders) {
      std::array<int, 4> ref {};
      for (int child = 0; child < 4; ++child)
        ref[child] = t8_mra::triangle_order::get_reference_children_order (type, child, order);
      EXPECT_TRUE (is_permutation<4> (ref))
        << "type " << type << " order {" << order[0] << order[1] << order[2] << "} not a bijection";
    }
}

}  // namespace

#endif  // T8_ENABLE_MRA
